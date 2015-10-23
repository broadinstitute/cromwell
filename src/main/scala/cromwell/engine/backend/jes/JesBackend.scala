package cromwell.engine.backend.jes

import java.math.BigInteger
import java.nio.file.{Path, Paths}

import com.google.api.services.genomics.model.Parameter
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding._
import cromwell.binding.expression.{NoFunctions, WdlStandardLibraryFunctions}
import cromwell.binding.types.{WdlArrayType, WdlFileType, WdlType}
import cromwell.binding.values._
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend._
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.jes.authentication._
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.engine.{AbortFunction, AbortRegistrationFunction, WorkflowDescriptor}
import cromwell.parser.BackendType
import cromwell.util.google.GoogleCloudStoragePath
import cromwell.util.{DigestionUtil, TryUtil}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object JesBackend {
  /*
    FIXME: At least for now the only files that can be used are stdout/stderr. However this leads to a problem
    where stdout.txt is input and output. Redirect stdout/stderr to a different name, but it'll be localized back
    in GCS as stdout/stderr. Yes, it's hacky.
   */
  val LocalWorkingDiskValue = s"disk://${RuntimeAttributes.LocalDiskName}"
  val ExecParamName = "exec"
  val WorkingDiskParamName = "working_disk"
  val ExtraConfigParamName = "__extra_config_gcs_path"
  val JesCromwellRoot = "/cromwell_root"
  val JesExecScript = "exec.sh"

  // Workflow options keys
  val RefreshTokenOptionKey = "refresh_token"
  val GcsRootOptionKey = "jes_gcs_root"
  val OptionKeys = Set(RefreshTokenOptionKey, GcsRootOptionKey)

  def authGcsCredentialsPath(gcsPath: String): JesInput = JesInput(ExtraConfigParamName, gcsPath, Paths.get(""), "LITERAL")

  // Decoration around WorkflowDescriptor to generate bucket names and the like
  implicit class JesWorkflowDescriptor(val descriptor: WorkflowDescriptor) extends JesBackend {
    def callDir(key: CallKey) = callGcsPath(descriptor, key.scope.name, key.index)
  }

  /**
   * Takes a path in GCS and comes up with a local path which is unique for the given GCS path
   * @param gcsPath The input path
   * @return A path which is unique per input path
   */
  def localFilePathFromCloudStoragePath(gcsPath: GoogleCloudStoragePath): Path = {
    Paths.get(JesCromwellRoot + "/" + gcsPath.bucket + "/" + gcsPath.objectName)
  }

  /**
   * Takes a possibly relative path and returns an absolute path, possibly under the JesCromwellRoot.
   * @param path The input path
   * @return A path which absolute
   */
  def localFilePathFromRelativePath(path: String): Path = {
    Paths.get(if (path.startsWith("/")) path else JesCromwellRoot + "/" + path)
  }

  /**
   * Takes a single WdlValue and maps google cloud storage (GCS) paths into an appropriate local file path.
   * If the input is not a WdlFile, or the WdlFile is not a GCS path, the mapping is a noop.
   *
   * @param wdlValue the value of the input
   * @return a new FQN to WdlValue pair, with WdlFile paths modified if appropriate.
   */
  private def gcsPathToLocal(wdlValue: WdlValue): WdlValue = {
    wdlValue match {
      case WdlFile(path, isGlob) =>
        GoogleCloudStoragePath.parse(path) match {
          case Success(gcsPath) => WdlFile(localFilePathFromCloudStoragePath(gcsPath).toString, isGlob)
          case Failure(e) => wdlValue
        }
      case array: WdlArray => array map gcsPathToLocal
      case _ => wdlValue
    }
  }

  protected def withRetry[T](f: Option[T] => T, failureMessage: String) = TryUtil.retryBlock(
    fn = f,
    retries = Some(10),
    pollingInterval = 5 seconds,
    pollingBackOffFactor = 1,
    maxPollingInterval = 10 seconds,
    failMessage = Some(failureMessage)
  )

  sealed trait JesParameter {
    def name: String
    def gcs: String
    def local: Path
    def parameterType: String

    final val toGoogleParameter = new Parameter().setName(name).setValue(local.toString).setType(parameterType)
  }

  final case class JesInput(name: String, gcs: String, local: Path, parameterType: String = "REFERENCE") extends JesParameter
  final case class JesOutput(name: String, gcs: String, local: Path, parameterType: String = "REFERENCE") extends JesParameter
}

class JesBackend extends Backend with LazyLogging with ProductionJesAuthentication with ProductionJesConfiguration {

  type BackendCall = JesBackendCall

  override def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = inputs map { case (k, v) => (k, gcsPathToLocal(v)) }
  override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs map { case (k,v) => (k, gcsPathToLocal(v)) }

  private def writeAuthenticationFile(workflow: WorkflowDescriptor) = authenticated { connection =>
    val path = GoogleCloudStoragePath(gcsAuthFilePath(workflow))

    generateAuthJson(jesConf.dockerCredentials, getGcsAuthInformation(workflow)) foreach { content =>
      def upload(prev: Option[Unit]) = connection.storage.uploadJson(path, content)

      logger.info(s"Creating authentication file for workflow ${workflow.id} at \n ${path.toString}")
      withRetry(upload, s"${makeTag(workflow)} Exception occurred while uploading auth file to $path")
    }
  }

  /**
   * Get a GcsLocalizing from workflow options if client secrets and refresh token are available.
   */
  def getGcsAuthInformation(workflow: WorkflowDescriptor): Option[JesAuthInformation] = {
    for {
      secrets <- jesConf.googleSecrets if jesConf.localizeWithRefreshToken
      token <- workflow.workflowOptions.get(RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(secrets, token)
  }

  /*
   * No need to copy GCS inputs for the workflow we should be able to directly reference them
   * Create an authentication json file containing docker credentials and/or user account information
   */
  override def initializeForWorkflow(workflow: WorkflowDescriptor): Try[HostInputs] = {
    writeAuthenticationFile(workflow)
    Success(workflow.actualInputs)
  }

  override def assertWorkflowOptions(options: WorkflowOptions): Unit = {
    // Warn for unrecognized option keys
    options.toMap.keySet.diff(OptionKeys) match {
      case unknowns if unknowns.nonEmpty => logger.warn(s"Unrecognized workflow option(s): ${unknowns.mkString(", ")}")
      case _ =>
    }

    if (jesConf.localizeWithRefreshToken) {
      Seq(RefreshTokenOptionKey) filterNot options.toMap.keySet match {
        case missing if missing.nonEmpty =>
          throw new IllegalArgumentException(s"Missing parameters in workflow options: ${missing.mkString(", ")}")
        case _ =>
      }
    }
  }

  /**
   * Delete authentication file in GCS once workflow is in a terminal state.
   *
   * First queries for the existence of the auth file, then deletes it if it exists.
   * If either of these operations fails, then a Future.failure is returned
   */
  override def cleanUpForWorkflow(workflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit] = authenticated { connection =>
    val gcsAuthFile = GoogleCloudStoragePath(gcsAuthFilePath(workflow))
    def gcsCheckAuthFileExists(prior: Option[Boolean]): Boolean = connection.storage.exists(gcsAuthFile)
    def gcsAttemptToDeleteObject(prior: Option[Unit]): Unit = connection.storage.deleteObject(gcsAuthFile)
    val tag = s"Workflow ${workflow.shortId}:"
    withRetry(gcsCheckAuthFileExists, s"$tag failed to query for auth file: $gcsAuthFile") match {
      case Success(exists) if exists =>
        withRetry(gcsAttemptToDeleteObject, s"$tag failed to delete auth file: $gcsAuthFile") match {
          case Success(_) => Future.successful(Unit)
          case Failure(ex) =>
            logger.error(s"$tag Could not delete the auth file $gcsAuthFile", ex)
            Future.failed(ex)
        }
      case Failure(ex) =>
        logger.error(s"$tag Could not query for the existence of the auth file $gcsAuthFile", ex)
        Future.failed(ex)
      case _ => Future.successful(Unit)
    }
  }

  override def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): StdoutStderr = {
    JesBackendCall.stdoutStderr(callGcsPath(descriptor, callName, index))
  }

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        key: CallKey,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortRegistrationFunction): BackendCall = {
    new JesBackendCall(this, workflowDescriptor, key, locallyQualifiedInputs, abortRegistrationFunction)
  }

  override def engineFunctions: WdlStandardLibraryFunctions = new JesEngineFunctionsWithoutCallContext()

  def execute(backendCall: BackendCall): ExecutionResult = {
    val tag = makeTag(backendCall)
    logger.info(s"$tag Call GCS path: ${backendCall.callGcsPath}")

    val jesInputs: Seq[JesInput] = generateJesInputs(backendCall).toSeq :+ backendCall.cmdInput
    val jesOutputs: Seq[JesOutput] = generateJesOutputs(backendCall)

    backendCall.instantiateCommand match {
      case Success(command) => runWithJes(backendCall, command, jesInputs, jesOutputs)
      case Failure(ex) => FailedExecution(ex)
    }
  }

  /**
   * Creates a set of JES inputs for a backend call.
   */
  def generateJesInputs(backendCall: BackendCall): Iterable[JesInput] = {
    val adjustedPaths = adjustInputPaths(backendCall.call, backendCall.locallyQualifiedInputs)
    def mkInput(locallyQualifiedInputName: String, location: WdlFile) = JesInput(locallyQualifiedInputName, backendCall.lookupFunction(locallyQualifiedInputName).valueString, Paths.get(location.valueString))
    adjustedPaths collect {
      case (locallyQualifiedInputName: String, location: WdlFile) =>
        logger.info(s"${makeTag(backendCall)}: $locallyQualifiedInputName -> ${location.valueString}")
        Seq(mkInput(locallyQualifiedInputName, location))
      case (locallyQualifiedInputName: String, localPathWdlArray: WdlArray) =>
        backendCall.lookupFunction(locallyQualifiedInputName) match {
          case WdlArray(WdlArrayType(memberType), remotePathArray: Seq[WdlValue]) if memberType == WdlFileType || memberType == WdlArrayType =>
            jesInputsFromArray(locallyQualifiedInputName, remotePathArray, localPathWdlArray.value)
          case _ => Seq[JesInput]() // Empty list.
        }
      case (locallyQualifiedInputName: String, location: WdlMap) => location.value flatMap { case (k, v) => Seq(k, v) } collect { case f: WdlFile => mkInput(locallyQualifiedInputName, f) }
    } flatten
  }

  /**
   * Takes two arrays of WDL values and generates any necessary JES inputs from them.
   */
  private def jesInputsFromArray(jesNamePrefix: String, remotePathArray: Seq[WdlValue], localPathArray: Seq[WdlValue]): Iterable[JesInput] = {
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((WdlArray(_, innerRemotePathArray: Seq[WdlValue]), WdlArray(_, innerLocalPathArray: Seq[WdlValue])), index) =>
        jesInputsFromArray(s"$jesNamePrefix-$index", innerRemotePathArray, innerLocalPathArray)
      case ((remotePath, localPath), index) => Seq(JesInput(s"$jesNamePrefix-$index", remotePath.valueString, Paths.get(localPath.valueString)))
    }
  }

  def generateJesOutputs(backendCall: BackendCall): Seq[JesOutput] = {
    val wdlFileOutputs = backendCall.call.task.outputs flatMap { taskOutput =>
      taskOutput.expression.evaluateFiles(backendCall.lookupFunction, new NoFunctions, taskOutput.wdlType) match {
        case Success(wdlFiles) => wdlFiles map gcsPathToLocal
        case Failure(ex) =>
          logger.warn(s"${makeTag(backendCall)} Could not evaluate $taskOutput: ${ex.getMessage}")
          Seq.empty[String]
      }
    }

    // Create the mappings. GLOB mappings require special treatment (i.e. stick everything matching the glob in a folder)
    wdlFileOutputs.distinct map {
      file => file match {
        case WdlFile(filePath, GlobOrNot.NotGlob) => JesOutput(makeSafeJesReferenceName(filePath), s"${backendCall.callGcsPath}/$filePath", localFilePathFromRelativePath(filePath))
        case WdlFile(filePath, GlobOrNot.YesGlob) => JesOutput(makeSafeJesReferenceName(filePath), s"${backendCall.callGcsPath}/globbed-${DigestionUtil.md5Sum(filePath)}/", localFilePathFromRelativePath(filePath))
      }
    }
  }

  /**
   * If the desired reference name is too long, we don't want to break JES or risk collisions by arbitrary truncation. So,
   * just use a hash. We only do this when needed to give better traceability in the normal case.
   */
  private def makeSafeJesReferenceName(referenceName: String) = {
    if (referenceName.length <= 127) referenceName else DigestionUtil.md5Sum(referenceName)
  }

  private def uploadCommandScript(backendCall: BackendCall, command: String): Try[Unit] = authenticated { implicit connection =>
    val fileContent =
      s"""
         |#!/bin/bash
         |cd $JesCromwellRoot && \\
         |$command
         |echo $$? > ${JesBackendCall.RcFilename}
       """.stripMargin.trim

    def attemptToUploadObject(priorAttempt: Option[Unit]) = authenticated { _.storage.uploadObject(backendCall.gcsExecPath, fileContent) }

    withRetry(attemptToUploadObject, s"${makeTag(backendCall)} Exception occurred while uploading script to ${backendCall.gcsExecPath}")
  }

  private def createJesRun(backendCall: BackendCall, jesParameters: Seq[JesParameter]): Try[Run] = authenticated { connection =>
    def attemptToCreateJesRun(priorAttempt: Option[Run]): Run = Pipeline(
      backendCall.jesCommandLine,
      backendCall.workflowDescriptor,
      backendCall.key,
      jesParameters,
      googleProject(backendCall.workflowDescriptor),
      connection
    ).run

    withRetry(attemptToCreateJesRun, s"${makeTag(backendCall)} Exception occurred while creating JES Run")
  }

  /**
   * Turns a GCS path representing a workflow input into the GCS path where the file would be mirrored to in this workflow:
   * task x {
   *  File x
   *  ...
   *  Output {
   *    File mirror = x
   *  }
   * }
   *
   * This function is more useful in working out the common prefix when the filename is modified somehow
   * in the workflow (e.g. "-new.txt" is appended)
   */
  private def gcsInputToGcsOutput(backendCall: BackendCall, inputValue: WdlValue): WdlValue = {
    // Convert to the local path where the file is localized to in the VM:
    val vmLocalizationPath = gcsPathToLocal(inputValue)

    vmLocalizationPath match {
      // If it's a file, work out where the file would be delocalized to, otherwise no-op:
      case x : WdlFile =>
        val delocalizationPath = s"${backendCall.callGcsPath}/${vmLocalizationPath.valueString}"
        WdlFile(delocalizationPath)
      case other => other
    }
  }

  private def customLookupFunction(backendCall: BackendCall) = { toBeLookedUp: String =>
    val originalLookup = backendCall.lookupFunction
    gcsInputToGcsOutput(backendCall, originalLookup(toBeLookedUp))
  }

  private def pollJesRun(run: Run, backendCall: BackendCall, jesOutputs: Seq[JesOutput]): ExecutionResult = authenticated { connection =>
    // Wait until the job starts (or completes/fails) before registering the abort to avoid awkward cancel-during-initialization behavior.
    val initializedStatus = run.waitUntilRunningOrComplete
    backendCall.callAbortRegistrationFunction.register(AbortFunction(() => run.abort()))

    def wdlFileToGcsPath(value: WdlValue, coerceTo: WdlType) =
      if (coerceTo == WdlFileType && WdlFileType.isCoerceableFrom(value.wdlType))
        jesOutputs find { _.name == value.valueString } map { j => WdlFile(j.gcs) } getOrElse value
      else
        value

    try {
      val status = run.waitUntilComplete(initializedStatus)

      val outputMappings = backendCall.call.task.outputs map { taskOutput =>
        /**
         * this will evaluate the task output expression and coerces it to the task output's type.
         * If the result is a WdlFile, then attempt to find the JesOutput with the same path and
         * return a WdlFile that represents the GCS path and not the local path.  For example,
         *
         * output {
         *   File x = "out" + ".txt"
         * }
         *
         * "out" + ".txt" is evaluated to WdlString("out.txt") and then coerced into a WdlFile("out.txt")
         * Then, via wdlFileToGcsPath(), we attempt to find the JesOutput with .name == "out.txt".
         * If it is found, then WdlFile("gs://some_bucket/out.txt") will be returned.
         */
        val attemptedValue = taskOutput.expression.evaluate(customLookupFunction(backendCall), backendCall.engineFunctions) map { wdlValue =>
          wdlFileToGcsPath(wdlValue, taskOutput.wdlType)
        }
        taskOutput.name -> attemptedValue
      } toMap

      lazy val stderrLength: BigInteger = connection.storage.objectSize(GoogleCloudStoragePath(backendCall.stderrJesOutput.gcs))
      lazy val returnCode = backendCall.downloadRcFile.map(_.trim.toInt)
      lazy val continueOnReturnCode = backendCall.call.continueOnReturnCode

      status match {
        case Run.Success if backendCall.call.failOnStderr && stderrLength.intValue > 0 =>
          FailedExecution(new Throwable(s"${makeTag(backendCall)} execution failed: stderr has length $stderrLength"))
        case Run.Success if returnCode.isFailure =>
          FailedExecution(new Throwable(s"${makeTag(backendCall)} execution failed: could not download or parse return code file", returnCode.failed.get))
        case Run.Success if !continueOnReturnCode.continueFor(returnCode.get) =>
          FailedExecution(new Throwable(s"${makeTag(backendCall)} execution failed: disallowed command return code: " + returnCode.get))
        case Run.Success =>
          handleSuccess(outputMappings, backendCall.workflowDescriptor, returnCode.get)
        case Run.Failed(errorCode, errorMessage) =>
          val throwable = if (errorMessage contains "Operation canceled at") {
            new TaskAbortedException()
          } else {
            new Throwable(s"Task ${backendCall.workflowDescriptor.id}:${backendCall.call.name} failed: error code $errorCode. Message: $errorMessage")
          }
          FailedExecution(throwable, Option(errorCode))

      }
    }
    catch {
      case e: Exception => FailedExecution(e)
    }
  }

  private def runWithJes(backendCall: BackendCall, command: String, jesInputs: Seq[JesInput], jesOutputs: Seq[JesOutput]): ExecutionResult = {
    val tag = makeTag(backendCall)
    val jesParameters = backendCall.standardParameters ++ gcsAuthParameter(backendCall.workflowDescriptor) ++ jesInputs ++ jesOutputs
    logger.info(s"$tag `$command`")

    val jesJobSetup = for {
      _ <- uploadCommandScript(backendCall, command)
      run <- createJesRun(backendCall, jesParameters)
    } yield run

    jesJobSetup match {
      case Failure(ex) =>
        logger.error(s"$tag Failed to create a JES run", ex)
        FailedExecution(ex)
      case Success(run) => pollJesRun(run, backendCall, jesOutputs)
    }
  }

  private def handleSuccess(outputMappings: Map[String, Try[WdlValue]], workflowDescriptor: WorkflowDescriptor, returnCode: Int): ExecutionResult = {
    val taskOutputEvaluationFailures = outputMappings filter { _._2.isFailure }
    val outputValues = if (taskOutputEvaluationFailures.isEmpty) {
      Success(outputMappings collect { case (name, Success(wdlValue)) => name -> wdlValue })
    } else {
      val message = taskOutputEvaluationFailures collect { case (name, Failure(e)) => s"$name: $e" } mkString "\n"
      Failure(new Throwable(s"Workflow ${workflowDescriptor.id}: $message"))
    }

    outputValues match {
      case Success(outputs) => SuccessfulExecution(outputs, returnCode)
      case Failure(e) => FailedExecution(e)
    }
  }

  override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow])(implicit ec: ExecutionContext): Future[Any] = Future("FIXME")

  override def backendType = BackendType.JES

  def gcsAuthFilePath(descriptor: WorkflowDescriptor): String = {
    val wfPath = workflowGcsPath(descriptor)
    s"$wfPath/auth.json"
  }

  def callGcsPath(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): String = {
    val shardPath = index map { i => s"/shard-$i"} getOrElse ""
    val workflowPath = workflowGcsPath(descriptor)
    s"$workflowPath/call-$callName$shardPath"
  }

  def workflowGcsPath(descriptor: WorkflowDescriptor): String = {
    val bucket = descriptor.workflowOptions.getOrElse(GcsRootOptionKey, jesConf.executionBucket)
    s"$bucket/${descriptor.namespace.workflow.name}/${descriptor.id}"
  }

  def googleProject(descriptor: WorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse("google_project", jesConf.project)
  }

  // Create an input parameter containing the path to this authentication file, if needed
  def gcsAuthParameter(descriptor: WorkflowDescriptor): Option[JesInput] = {
    if (jesConf.localizeWithRefreshToken || jesConf.isDockerAuthenticated)
      Option(authGcsCredentialsPath(gcsAuthFilePath(descriptor)))
    else None
  }
}
