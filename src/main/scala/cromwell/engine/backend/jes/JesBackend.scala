package cromwell.engine.backend.jes

import java.io.File
import java.math.BigInteger
import java.nio.file.{Path, Paths}

import com.google.api.services.genomics.model.Parameter
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding._
import cromwell.binding.expression.NoFunctions
import cromwell.binding.types.{WdlFileType, WdlType}
import cromwell.binding.values._
import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend._
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.engine.{AbortFunction, AbortRegistrationFunction, WorkflowDescriptor}
import cromwell.parser.BackendType
import cromwell.util.TryUtil
import cromwell.util.google.GoogleCloudStoragePath

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
  val LocalWorkingDiskValue = "disk://local-disk"
  val WorkingDiskParamName = "working_disk"
  val ExtraConfigParamName = "__extra_config_gcs_path"
  val JesCromwellRoot = "/cromwell_root"

  // Workflow options keys
  val AccountOptionKey = "account_name"
  val RefreshTokenOptionKey = "refresh_token"
  val GcsRootOptionKey = "jes_gcs_root"
  val OptionKeys = Set(AccountOptionKey, RefreshTokenOptionKey, GcsRootOptionKey)

  def localizationDiskInput(): JesInput = JesInput(WorkingDiskParamName, LocalWorkingDiskValue, new File(JesCromwellRoot).toPath)

  def authGcsCredentialsPath(gcsPath: Option[String]): Option[JesInput] =
    gcsPath.map(JesInput(ExtraConfigParamName, _, Paths.get(""), "LITERAL"))

  // Decoration around WorkflowDescriptor to generate bucket names and the like
  implicit class JesWorkflowDescriptor(val descriptor: WorkflowDescriptor) extends JesBackend {
    def callDir(key: CallKey) = callGcsPath(descriptor, key.scope.name, key.index)
  }

  sealed trait JesParameter {
    def name: String
    def gcs: String
    def local: Path
    def parameterType: String

    final val isInput = this.isInstanceOf[JesInput]
    final val isOutput = !isInput
    final val toGoogleParameter = new Parameter().setName(name).setValue(local.toString).setType(parameterType)
  }

  final case class JesInput(name: String, gcs: String, local: Path, parameterType: String = "REFERENCE") extends JesParameter
  final case class JesOutput(name: String, gcs: String, local: Path, parameterType: String = "REFERENCE") extends JesParameter
}

class JesBackend extends Backend with LazyLogging {
  type BackendCall = JesBackendCall

  /*
   * NOTE: Those have been moved from JesBackend object to the class,
   * so if we were ever to instantiate several backends, as many JesConf/JesConnection would be created as well
   */
  lazy val JesConf = JesAttributes()
  lazy val JesConnection = JesInterface(JesConf.applicationName, JesConf.endpointUrl)

  /**
   * Takes a path in GCS and comes up with a local path which is unique for the given GCS path
   * @param gcsPath The input path
   * @return A path which is unique per input path
   */
  def localFilePathFromCloudStoragePath(gcsPath: GoogleCloudStoragePath): Path = {
    Paths.get(JesCromwellRoot + "/" + gcsPath.bucket + "/" + gcsPath.objectName)
  }

  /**
   * Takes a single FQN to WdlValue pair and maps google cloud storage (GCS) paths into an appropriate local file path.
   * If the input is not a WdlFile, or the WdlFile is not a GCS path, the mapping is a noop.
   *
   * @param fqn the FQN of the input variable
   * @param wdlValue the value of the input
   * @return a new FQN to WdlValue pair, with WdlFile paths modified if appropriate.
   */
  def mapInputValue(fqn: String, wdlValue: WdlValue): (String, WdlValue) = {
    wdlValue match {
      case WdlFile(path) =>
        GoogleCloudStoragePath.parse(path) match {
          case Success(gcsPath) => (fqn, WdlFile(localFilePathFromCloudStoragePath(gcsPath).toString))
          case Failure(e) => (fqn, wdlValue)
        }
      case array:WdlArray => fqn -> array.map(v => mapInputValue("bogus", v)._2) // TODO: this is ugly.
      case _ => (fqn, wdlValue)
    }
  }

  override def adjustInputPaths(call: Call, inputs: CallInputs): CallInputs = inputs map { case (k, v) => mapInputValue(k, v) }
  override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs

  private def writeAuthenticationFile(workflow: WorkflowDescriptor, userAuth: Option[GcsUserAuthInformation]) = {
    val path = GoogleCloudStoragePath(gcsAuthFilePath(workflow))
    GcsAuth.generateJson(JesConf.dockerCredentials, userAuth) foreach { content =>
      logger.info(s"Creating authentication file for workflow ${workflow.id} at \n ${path.toString}")
      JesConnection.storage.uploadJson(path, content)
    }
  }

  /*
   * No need to copy GCS inputs for the workflow we should be able to directly reference them
   * Create an authentication json file containing docker credentials and/or user account information
   */
  override def initializeForWorkflow(workflow: WorkflowDescriptor): Try[HostInputs] = {
    def writeWithRefreshToken() = getGcsAuthInformation(workflow) map { userAuth => writeAuthenticationFile(workflow, Option(userAuth)) }

    JesConf.authMode match {
      case RefreshTokenMode => writeWithRefreshToken().map(_ => workflow.actualInputs)
      case _ =>
        writeAuthenticationFile(workflow, None)
        Success(workflow.actualInputs)
    }
  }

  override def assertWorkflowOptions(options: WorkflowOptions): Unit = {
    // Warn for unrecognized option keys
    options.toMap.keySet.diff(OptionKeys) match {
      case unknowns if unknowns.nonEmpty => logger.warn(s"Unrecognized workflow option(s): ${unknowns.mkString(", ")}")
      case _ =>
    }

    if (JesConf.authMode == RefreshTokenMode) {
      Seq(AccountOptionKey, RefreshTokenOptionKey) filterNot options.toMap.keySet match {
        case missing if missing.nonEmpty =>
          throw new IllegalArgumentException(s"Missing parameters in workflow options: ${missing.mkString(", ")}")
        case _ =>
      }
    }
  }

  /**
   * Delete authentication file in gcs once workflow is in a terminal state.
   */
  override def cleanUpForWorkflow(workflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit] = {
    Future(JesConnection.storage.deleteObject(GoogleCloudStoragePath(gcsAuthFilePath(workflow))))
  }

  /**
   * Get a GcsUserAuthInformation from workflow options
   */
  def getGcsAuthInformation(workflow: WorkflowDescriptor): Try[GcsUserAuthInformation] = {
    for {
      account <- workflow.workflowOptions.get(AccountOptionKey)
      token <- workflow.workflowOptions.get(RefreshTokenOptionKey)
    } yield GcsUserAuthInformation(account, token)
  }

  override def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): StdoutStderr = {
    JesBackendCall.stdoutStderr(callGcsPath(descriptor, callName, index))
  }

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        key: CallKey,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortRegistrationFunction): BackendCall = {
    JesBackendCall(this, workflowDescriptor, key, locallyQualifiedInputs, abortRegistrationFunction)
  }

  def execute(backendCall: BackendCall): ExecutionResult = {
    val tag = makeTag(backendCall)
    val cmdInput = JesInput("exec", backendCall.gcsExecPath.toString, Paths.get("exec.sh"))
    logger.info(s"$tag Call GCS path: ${backendCall.callGcsPath}")

    val jesInputs: Seq[JesInput] = generateJesInputs(backendCall).toSeq :+ cmdInput
    val jesOutputs: Seq[JesOutput] = generateJesOutputs(backendCall)

    backendCall.instantiateCommand match {
      case Success(command) => runWithJes(backendCall, command, jesInputs, jesOutputs)
      case Failure(ex) => FailedExecution(ex)
    }
  }

  def generateJesInputs(backendCall: BackendCall): Iterable[JesInput] = {
    val adjustedPaths = adjustInputPaths(backendCall.call, backendCall.locallyQualifiedInputs)
    def mkInput(fileTag: String, location: WdlFile) = JesInput(fileTag, backendCall.lookupFunction(fileTag).valueString, Paths.get(location.valueString))
    adjustedPaths collect {
      case (fileTag: String, location: WdlFile) =>
        logger.info(s"${makeTag(backendCall)}: $fileTag -> ${location.valueString}")
        Seq(mkInput(fileTag, location))
      case (fileTag: String, location: WdlArray) => location.value.collect { case f: WdlFile => mkInput(fileTag, f) }
      case (fileTag: String, location: WdlMap) => location.value flatMap { case (k, v) => Seq(k, v) } collect { case f: WdlFile => mkInput(fileTag, f) }
    } flatten
  }

  def generateJesOutputs(backendCall: BackendCall): Seq[JesOutput] = {
    val wdlFileOutputs = backendCall.call.task.outputs flatMap { taskOutput =>
      taskOutput.expression.evaluateFiles(backendCall.lookupFunction, new NoFunctions, taskOutput.wdlType) match {
        case Success(wdlFiles) => wdlFiles.map(_.valueString)
        case Failure(ex) =>
          logger.warn(s"${makeTag(backendCall)} Could not evaluate $taskOutput: ${ex.getMessage}")
          Seq.empty[String]
      }
    }
    wdlFileOutputs.distinct map { filePath =>
      JesOutput(filePath, s"${backendCall.callGcsPath}/$filePath", Paths.get(filePath))
    }
  }

  private def uploadCommandScript(backendCall: BackendCall, command: String): Try[Unit] = {
    val fileContent =
      s"""
         |#!/bin/bash
         |$command
         |echo $$? > ${JesBackendCall.RcFilename}
       """.stripMargin.trim

    def attemptToUploadObject(priorAttempt: Option[Unit]) = JesConnection.storage.uploadObject(backendCall.gcsExecPath, fileContent)

    TryUtil.retryBlock(
      fn = attemptToUploadObject,
      retries = Some(10),
      pollingInterval = 5 seconds,
      pollingBackOffFactor = 1,
      maxPollingInterval = 10 seconds,
      failMessage = Some(s"${makeTag(backendCall)} Exception occurred while uploading script to ${backendCall.gcsExecPath}")
    )
  }

  private def createJesRun(backendCall: BackendCall, jesParameters: Seq[JesParameter]): Try[Run] = {
    def attemptToCreateJesRun(priorAttempt: Option[Run]): Run = Pipeline(
      backendCall.jesCommandLine,
      backendCall.workflowDescriptor,
      backendCall.key,
      jesParameters,
      googleProject(backendCall.workflowDescriptor),
      JesConnection
    ).run

    TryUtil.retryBlock(
      fn = attemptToCreateJesRun,
      retries = Some(10),
      pollingInterval = 5 seconds,
      pollingBackOffFactor = 1,
      maxPollingInterval = 10 seconds,
      failMessage = Some(s"${makeTag(backendCall)} Exception occurred while creating JES Run")
    )
  }

  private def pollJesRun(run: Run, backendCall: BackendCall, jesOutputs: Seq[JesOutput]): ExecutionResult = {
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

      val outputMappings = backendCall.call.task.outputs map {taskOutput =>
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
        val attemptedValue = taskOutput.expression.evaluate(backendCall.lookupFunction, backendCall.engineFunctions) map { wdlValue =>
          wdlFileToGcsPath(wdlValue, taskOutput.wdlType)
        }
        taskOutput.name -> attemptedValue
      } toMap

      lazy val stderrLength: BigInteger = JesConnection.storage.objectSize(GoogleCloudStoragePath(backendCall.stderrJesOutput.gcs))
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
          FailedExecution(throwable)
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
    s"$wfPath/gcloudauth.json"
  }

  def callGcsPath(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): String = {
    val shardPath = index map { i => s"/shard-$i"} getOrElse ""
    val workflowPath = workflowGcsPath(descriptor)
    s"$workflowPath/call-$callName$shardPath"
  }

  def workflowGcsPath(descriptor: WorkflowDescriptor): String = {
    val bucket = descriptor.workflowOptions.getOrElse(GcsRootOptionKey, JesConf.executionBucket)
    s"$bucket/${descriptor.namespace.workflow.name}/${descriptor.id}"
  }

  def googleProject(descriptor: WorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse("google_project", JesConf.project)
  }

  // Create an input parameter containing the path to this authentication file
  def gcsAuthParameter(descriptor: WorkflowDescriptor): Option[JesInput] = {
    if (JesConf.authMode == RefreshTokenMode || JesConf.dockerCredentials.isDefined) {
      authGcsCredentialsPath(Option(gcsAuthFilePath(descriptor)))
    } else None
  }
}
