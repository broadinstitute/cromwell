package cromwell.engine.backend.jes

import java.math.BigInteger
import java.net.SocketTimeoutException
import java.nio.file.{Path, Paths}

import com.google.api.services.genomics.model.Parameter
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding._
import cromwell.binding.expression.{NoFunctions, WdlStandardLibraryFunctions}
import cromwell.binding.values._
import cromwell.engine.ExecutionIndex.{ExecutionIndex, IndexEnhancedInt}
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend._
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.jes.Run.RunStatus
import cromwell.engine.backend.jes.authentication._
import cromwell.engine.db.DataAccess.globalDataAccess
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.db.slick.Execution
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor, _}
import cromwell.logging.WorkflowLogger
import cromwell.parser.BackendType
import cromwell.util.StringUtil._
import cromwell.util.TryUtil
import cromwell.util.google.{GoogleCloudStorage, GoogleCloudStoragePath, GoogleCredentialFactory}

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
  val MonitoringParamName = "monitoring"
  val WorkingDiskParamName = "working_disk"
  val ExtraConfigParamName = "__extra_config_gcs_path"
  val JesCromwellRoot = "/cromwell_root"
  val JesExecScript = "exec.sh"
  val JesMonitoringScript = "monitoring.sh"
  val JesMonitoringLogFile = "monitoring.log"

  // Workflow options keys
  val RefreshTokenOptionKey = "refresh_token"
  val GcsRootOptionKey = "jes_gcs_root"
  val MonitoringScriptOptionKey = "monitoring_script"
  val GoogleProjectOptionKey = "google_project"
  val AuthFilePathOptionKey = "auth_bucket"
  val OptionKeys = Set(RefreshTokenOptionKey, GcsRootOptionKey, MonitoringScriptOptionKey, GoogleProjectOptionKey, AuthFilePathOptionKey)

  lazy val monitoringScriptLocalPath = localFilePathFromRelativePath(JesMonitoringScript)
  lazy val monitoringLogFileLocalPath = localFilePathFromRelativePath(JesMonitoringLogFile)

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
      case wdlFile: WdlFile =>
        GoogleCloudStoragePath.parse(wdlFile.value) match {
          case Success(gcsPath) => WdlFile(localFilePathFromCloudStoragePath(gcsPath).toString, wdlFile.isGlob)
          case Failure(e) => wdlValue
        }
      case wdlArray: WdlArray => wdlArray map gcsPathToLocal
      case wdlMap: WdlMap => wdlMap map { case (k, v) => gcsPathToLocal(k) -> gcsPathToLocal(v) }
      case _ => wdlValue
    }
  }

  protected def withRetry[T](f: Option[T] => T, logger: WorkflowLogger, failureMessage: String) = TryUtil.retryBlock(
    fn = f,
    retryLimit = Option(10),
    pollingInterval = 5 seconds,
    pollingBackOffFactor = 1,
    maxPollingInterval = 10 seconds,
    logger = logger,
    failMessage = Option(failureMessage)
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

  implicit class EnhancedExecution(val execution: Execution) extends AnyVal {
    import cromwell.engine.ExecutionIndex._
    def toKey: ExecutionDatabaseKey = ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex)
    def isScatter: Boolean = execution.callFqn.contains(Scatter.FQNIdentifier)
    def executionStatus: ExecutionStatus = ExecutionStatus.withName(execution.status)
  }

  type IOInterface = GoogleCloudStorage
}

final case class JesJobKey(operationId: String) extends JobKey

/**
 * Representing a running JES execution, instances of this class are never Done and it is never okay to
 * ask them for results.
 */
case class JesPendingExecutionHandle(backendCall: JesBackendCall,
                                     jesOutputs: Seq[JesOutput],
                                     run: Run,
                                     previousStatus: Option[RunStatus]) extends ExecutionHandle {
  override val isDone = false

  override def result = FailedExecution(new IllegalStateException)
}


class JesBackend extends Backend with LazyLogging with ProductionJesAuthentication with ProductionJesConfiguration {

  type BackendCall = JesBackendCall
  type IOInterface = JesBackend.IOInterface

  override def adjustInputPaths(callKey: CallKey, inputs: CallInputs, workflowDescriptor: WorkflowDescriptor): CallInputs = inputs mapValues gcsPathToLocal
  override def adjustOutputPaths(call: Call, outputs: CallOutputs): CallOutputs = outputs mapValues {
    case CallOutput(value, hash) => CallOutput(gcsPathToLocal(value), hash) }

  def fileHasher(workflow: WorkflowDescriptor): FileHasher = { wdlFile: WdlFile => SymbolHash(getCrc32c(workflow, GoogleCloudStoragePath(wdlFile.value))) }

  private def writeAuthenticationFile(workflow: WorkflowDescriptor) = authenticateAsCromwell { connection =>
    val log = workflowLogger(workflow)

    generateAuthJson(jesConf.dockerCredentials, getGcsAuthInformation(workflow)) foreach { content =>
      val path = GoogleCloudStoragePath(gcsAuthFilePath(workflow))
      def upload(prev: Option[Unit]) = connection.storage.uploadJson(path, content)

      log.info(s"Creating authentication file for workflow ${workflow.id} at \n ${path.toString}")
      withRetry(upload, log, s"Exception occurred while uploading auth file to $path")
    }
  }

  def getCrc32c(workflow: WorkflowDescriptor, googleCloudStoragePath: GoogleCloudStoragePath): String = authenticateAsUser(workflow) {
    _.getCrc32c(googleCloudStoragePath)
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
  override def cleanUpForWorkflow(workflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit] = {
    Future(gcsAuthFilePath(workflow)) map { path =>
      deleteAuthFile(path, workflowLogger(workflow))
      ()
    } recover {
      case e: UnsatisfiedInputsException =>  // No need to fail here, it just means that we didn't have an auth file in the first place so no need to delete it.
    }
  }

  private def deleteAuthFile(authFilePath: String, log: WorkflowLogger): Future[Unit] = authenticateAsCromwell { connection =>
      def gcsCheckAuthFileExists(prior: Option[Boolean]): Boolean = connection.storage.exists(authFilePath)
      def gcsAttemptToDeleteObject(prior: Option[Unit]): Unit = connection.storage.deleteObject(authFilePath)
      withRetry(gcsCheckAuthFileExists, log, s"Failed to query for auth file: $authFilePath") match {
        case Success(exists) if exists =>
          withRetry(gcsAttemptToDeleteObject, log, s"Failed to delete auth file: $authFilePath") match {
            case Success(_) => Future.successful(Unit)
            case Failure(ex) =>
              log.error(s"Could not delete the auth file $authFilePath", ex)
              Future.failed(ex)
          }
        case Failure(ex) =>
          log.error(s"Could not query for the existence of the auth file $authFilePath", ex)
          Future.failed(ex)
        case _ => Future.successful(Unit)
      }
  }

  override def stdoutStderr(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): CallLogs = {
    JesBackendCall.stdoutStderr(callGcsPath(descriptor, callName, index))
  }

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        key: CallKey,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortRegistrationFunction): BackendCall = {
    new JesBackendCall(this, workflowDescriptor, key, locallyQualifiedInputs, abortRegistrationFunction)
  }

  override def engineFunctions(interface: IOInterface): WdlStandardLibraryFunctions = new JesEngineFunctionsWithoutCallContext(interface)

  private def executeOrResume(backendCall: BackendCall, runIdForResumption: Option[String])(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    val log = workflowLoggerWithCall(backendCall)
    log.info(s"Call GCS path: ${backendCall.callGcsPath}")
    val monitoringScript: Option[JesInput] = monitoringIO(backendCall)
    val monitoringOutput = monitoringScript map { _ => JesOutput(MonitoringParamName, backendCall.defaultMonitoringOutputPath, monitoringLogFileLocalPath) }

    val jesInputs: Seq[JesInput] = generateJesInputs(backendCall).toSeq ++ monitoringScript :+ backendCall.cmdInput
    val jesOutputs: Seq[JesOutput] = generateJesOutputs(backendCall) ++ monitoringOutput

    backendCall.instantiateCommand match {
      case Success(command) => runWithJes(backendCall, command, jesInputs, jesOutputs, runIdForResumption, monitoringScript.isDefined)
      case Failure(ex: SocketTimeoutException) => throw ex // probably a GCS transient issue, throwing will cause it to be retried
      case Failure(ex) => FailedExecutionHandle(ex)
    }
  }

  def execute(backendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = executeOrResume(backendCall, runIdForResumption = None)

  def resume(backendCall: BackendCall, jobKey: JobKey)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    val runId = Option(jobKey) collect { case jesKey: JesJobKey => jesKey.operationId }
    executeOrResume(backendCall, runIdForResumption = runId)
  }

  /**
   * Creates a set of JES inputs for a backend call.
   * Note that duplicates input files (same gcs path) will be (re-)localized every time they are referenced.
   */
  def generateJesInputs(backendCall: BackendCall): Iterable[JesInput] = {
    backendCall.locallyQualifiedInputs mapValues { _.collectAsSeq { case w: WdlFile => w } } flatMap {
      case (name, files) => jesInputsFromWdlFiles(name, files, files map { gcsPathToLocal(_).asInstanceOf[WdlFile] })
    }
  }

  /**
   * Takes two arrays of remote and local WDL File paths and generates the necessary JESInputs.
   */
  private def jesInputsFromWdlFiles(jesNamePrefix: String, remotePathArray: Seq[WdlFile], localPathArray: Seq[WdlFile]): Iterable[JesInput] = {
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((remotePath, localPath), index) => Seq(JesInput(s"$jesNamePrefix-$index", remotePath.valueString, Paths.get(localPath.valueString)))
    }
  }

  def generateJesOutputs(backendCall: BackendCall): Seq[JesOutput] = {
    val log = workflowLoggerWithCall(backendCall)
    val wdlFileOutputs = backendCall.call.task.outputs flatMap { taskOutput =>
      taskOutput.expression.evaluateFiles(backendCall.lookupFunction, new NoFunctions, taskOutput.wdlType) match {
        case Success(wdlFiles) => wdlFiles map gcsPathToLocal
        case Failure(ex) =>
          log.warn(s"Could not evaluate $taskOutput: ${ex.getMessage}")
          Seq.empty[String]
      }
    }

    // Create the mappings. GLOB mappings require special treatment (i.e. stick everything matching the glob in a folder)
    wdlFileOutputs.distinct map {
      case wdlFile: WdlFile =>
        val destination = wdlFile match {
          case WdlSingleFile(filePath) => s"${backendCall.callGcsPath}/$filePath"
          case WdlGlobFile(filePath) => backendCall.globOutputPath(filePath)
        }
        JesOutput(makeSafeJesReferenceName(wdlFile.value), destination, localFilePathFromRelativePath(wdlFile.value))
    }
  }

  def monitoringIO(backendCall: BackendCall): Option[JesInput] = {
    backendCall.workflowDescriptor.workflowOptions.get(MonitoringScriptOptionKey) map { path =>
      JesInput(MonitoringParamName, GoogleCloudStoragePath(path).toString, monitoringScriptLocalPath)
    } toOption
  }

  /**
   * If the desired reference name is too long, we don't want to break JES or risk collisions by arbitrary truncation. So,
   * just use a hash. We only do this when needed to give better traceability in the normal case.
   */
  private def makeSafeJesReferenceName(referenceName: String) = {
    if (referenceName.length <= 127) referenceName else referenceName.md5Sum
  }

  private def uploadCommandScript(backendCall: BackendCall, command: String, withMonitoring: Boolean): Try[Unit] = {
    val monitoring = if (withMonitoring) {
      s"""|touch $JesMonitoringLogFile
          |chmod u+x $JesMonitoringScript
          |./$JesMonitoringScript > $JesMonitoringLogFile &""".stripMargin
    } else ""

    val fileContent =
      s"""
         |#!/bin/bash
         |export _JAVA_OPTIONS=-Djava.io.tmpdir=$JesCromwellRoot/tmp
         |export TMPDIR=$JesCromwellRoot/tmp
         |cd $JesCromwellRoot
         |$monitoring
         |$command
         |echo $$? > ${JesBackendCall.RcFilename}
       """.stripMargin.trim

    def attemptToUploadObject(priorAttempt: Option[Unit]) = authenticateAsUser(backendCall.workflowDescriptor) { _.uploadObject(backendCall.gcsExecPath, fileContent) }

    val log = workflowLogger(backendCall.workflowDescriptor)
    withRetry(attemptToUploadObject, log, s"${workflowLoggerWithCall(backendCall).tag} Exception occurred while uploading script to ${backendCall.gcsExecPath}")
  }

  private def createJesRun(backendCall: BackendCall, jesParameters: Seq[JesParameter], runIdForResumption: Option[String]): Try[Run] =
    authenticateAsCromwell { connection =>
      def attemptToCreateJesRun(priorAttempt: Option[Run]): Run = Pipeline(
        backendCall.jesCommandLine,
        backendCall.workflowDescriptor,
        backendCall.key,
        jesParameters,
        googleProject(backendCall.workflowDescriptor),
        connection,
        runIdForResumption
      ).run

      val log = workflowLogger(backendCall.workflowDescriptor)
      withRetry(attemptToCreateJesRun, log, s"${workflowLoggerWithCall(backendCall).tag} Exception occurred while creating JES Run")
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
      case a: WdlArray => WdlArray(a.wdlType, a.value map { f => gcsInputToGcsOutput(backendCall, f) })
      case m: WdlMap => WdlMap(m.wdlType, m.value map { case (k, v) => gcsInputToGcsOutput(backendCall, k) -> gcsInputToGcsOutput(backendCall, v) })
      case other => other
    }
  }

  private def customLookupFunction(backendCall: BackendCall) = { toBeLookedUp: String =>
    val originalLookup = backendCall.lookupFunction
    gcsInputToGcsOutput(backendCall, originalLookup(toBeLookedUp))
  }

  def wdlValueToGcsPath(jesOutputs: Seq[JesOutput])(value: WdlValue): WdlValue = {
    def toGcsPath(wdlFile: WdlFile) = jesOutputs collectFirst { case o if o.name == makeSafeJesReferenceName(wdlFile.valueString) => WdlFile(o.gcs) } getOrElse value
    value match {
      case wdlArray: WdlArray => wdlArray map wdlValueToGcsPath(jesOutputs)
      case wdlMap: WdlMap => wdlMap map {
        case (k, v) => wdlValueToGcsPath(jesOutputs)(k) -> wdlValueToGcsPath(jesOutputs)(v)
      }
      case file: WdlFile => toGcsPath(file)
      case other => other
    }
  }

  def executionResult(status: RunStatus, handle: JesPendingExecutionHandle): ExecutionHandle = {
    val log = workflowLoggerWithCall(handle.backendCall)

    try {
      val backendCall = handle.backendCall
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
        val attemptedValue = taskOutput.expression.evaluate(customLookupFunction(backendCall), backendCall.engineFunctions) flatMap { wdlValue =>
          taskOutput.wdlType.coerceRawValue(wdlValue) map wdlValueToGcsPath(handle.jesOutputs)
        }
        taskOutput.name -> attemptedValue
      } toMap

      lazy val stderrLength: BigInteger = authenticateAsUser(backendCall.workflowDescriptor) { _.objectSize(GoogleCloudStoragePath(backendCall.stderrJesOutput.gcs)) }
      lazy val returnCodeContents = backendCall.downloadRcFile
      lazy val returnCode = returnCodeContents map { _.trim.toInt }
      lazy val continueOnReturnCode = backendCall.call.continueOnReturnCode

      status match {
        case Run.Success if backendCall.call.failOnStderr && stderrLength.intValue > 0 =>
          FailedExecutionHandle(new Throwable(s"${log.tag} execution failed: stderr has length $stderrLength"))
        case Run.Success if returnCodeContents.isFailure =>
          val exception = returnCode.failed.get
          log.warn(s"${log.tag} could not download return code file, retrying: " + exception.getMessage, exception)
          // Return handle to try again.
          handle
        case Run.Success if returnCode.isFailure =>
          FailedExecutionHandle(new Throwable(s"${log.tag} execution failed: could not parse return code as integer: " + returnCodeContents.get))
        case Run.Success if !continueOnReturnCode.continueFor(returnCode.get) =>
          FailedExecutionHandle(new Throwable(s"${log.tag} execution failed: disallowed command return code: " + returnCode.get))
        case Run.Success =>
          handleSuccess(outputMappings, backendCall.workflowDescriptor, returnCode.get, handle)
        case Run.Failed(errorCode, errorMessage) =>
          val throwable = if (errorMessage contains "Operation canceled at") {
            new TaskAbortedException()
          } else {
            new Throwable(s"Task ${backendCall.workflowDescriptor.id}:${backendCall.call.name} failed: error code $errorCode. Message: $errorMessage")
          }
          FailedExecutionHandle(throwable, Option(errorCode))
      }
    } catch {
      case e: Exception =>
        log.warn("Caught exception trying to download result, retrying: " + e.getMessage, e)
        // Return the original handle to try again.
        handle
    }
  }

  private def runWithJes(backendCall: BackendCall,
                         command: String,
                         jesInputs: Seq[JesInput],
                         jesOutputs: Seq[JesOutput],
                         runIdForResumption: Option[String],
                         withMonitoring: Boolean): ExecutionHandle = {
    val log = workflowLoggerWithCall(backendCall)
    val jesParameters = backendCall.standardParameters ++ gcsAuthParameter(backendCall.workflowDescriptor) ++ jesInputs ++ jesOutputs
    log.info(s"`$command`")

    val jesJobSetup = for {
      _ <- uploadCommandScript(backendCall, command, withMonitoring)
      run <- createJesRun(backendCall, jesParameters, runIdForResumption)
    } yield run

    jesJobSetup match {
      case Failure(ex) =>
        log.warn(s"Failed to create a JES run", ex)
        throw ex  // Probably a transient issue, throwing retries it
      case Success(run) => JesPendingExecutionHandle(backendCall, jesOutputs, run, previousStatus = None)
    }
  }

  private def handleSuccess(outputMappings: Map[String, Try[WdlValue]],
                            workflowDescriptor: WorkflowDescriptor,
                            returnCode: Int,
                            executionHandle: ExecutionHandle): ExecutionHandle = {

    val taskOutputEvaluationFailures = outputMappings filter { _._2.isFailure }
    if (taskOutputEvaluationFailures.isEmpty) {
      val outputs = outputMappings collect { case (name, Success(wdlValue)) => name -> CallOutput(wdlValue, wdlValue.getHash(fileHasher(workflowDescriptor))) }
      SuccessfulExecutionHandle(outputs, returnCode)
    } else if (taskOutputEvaluationFailures forall (_._2.failed.get.isInstanceOf[SocketTimeoutException])) {
      // Assume this as a transient exception trying to do some expression evaluation that involves reading from GCS.
      // Return the input handle to try again.
      executionHandle
    } else {
      val message = taskOutputEvaluationFailures collect { case (name, Failure(e)) => s"$name: $e" } mkString "\n"
      FailedExecutionHandle(new Throwable(s"Workflow ${workflowDescriptor.id}: $message"))
    }
  }

  /**
   * <ul>
   *   <li>Any execution in Failed should fail the restart.</li>
   *   <li>Any execution in Aborted should fail the restart.</li>
   *   <li>Scatters in Starting should fail the restart.</li>
   *   <li>Collectors in Running should be set back to NotStarted.</li>
   *   <li>Calls in Starting should be rolled back to NotStarted.</li>
   *   <li>Calls in Running with no job key should be rolled back to NotStarted.</li>
   * </ul>
   *
   * Calls in Running *with* a job key should be left in Running.  The WorkflowActor is responsible for
   * resuming the CallActors for these calls.
   */
  override def prepareForRestart(restartableWorkflow: WorkflowDescriptor)(implicit ec: ExecutionContext): Future[Unit] = {

    lazy val tag = s"Workflow ${restartableWorkflow.id.shortString}:"

    def handleExecutionStatuses(executions: Traversable[Execution]): Future[Unit] = {

      def stringifyExecutions(executions: Traversable[Execution]): String = {
        executions.toSeq.sortWith((lt, rt) => lt.callFqn < rt.callFqn || (lt.callFqn == rt.callFqn && lt.index < rt.index)).mkString(" ")
      }

      def isRunningCollector(key: Execution) = key.index.toIndex.isEmpty && key.executionStatus == ExecutionStatus.Running

      val failedOrAbortedExecutions = executions filter { x => x.executionStatus == ExecutionStatus.Aborted || x.executionStatus == ExecutionStatus.Failed }

      if (failedOrAbortedExecutions.nonEmpty) {
        Future.failed(new Throwable(s"$tag Cannot restart, found Failed and/or Aborted executions: " + stringifyExecutions(failedOrAbortedExecutions)))
      } else {
        // Cromwell has execution types: scatter, collector, call.
        val (scatters, collectorsAndCalls) = executions partition { _.isScatter }
        // If a scatter is found in starting state, it's not clear without further database queries whether the call
        // shards have been created or not.  This is an unlikely scenario and could be worked around with further
        // queries or a bracketing transaction, but for now Cromwell just bails out on restarting the workflow.
        val startingScatters = scatters filter { _.executionStatus == ExecutionStatus.Starting }
        if (startingScatters.nonEmpty) {
          Future.failed(new Throwable(s"$tag Cannot restart, found scatters in Starting status: " + stringifyExecutions(startingScatters)))
        } else {
          // Scattered calls have more than one execution with the same FQN.  Find any collectors in these FQN
          // groupings which are in Running state.
          // This is a race condition similar to the "starting scatters" case above, but here the assumption is that
          // it's more likely that collectors can safely be reset to starting.  This may prove not to be the case if
          // entries have been written to the symbol table.
          // Like the starting scatters case, further queries or a bracketing transaction would be a better long term solution.
          val runningCollectors = collectorsAndCalls.groupBy(_.callFqn) collect {
            case (_, xs) if xs.size > 1 => xs filter isRunningCollector } flatten

          for {
            _ <- globalDataAccess.resetNonResumableJesExecutions(restartableWorkflow.id)
            _ <- globalDataAccess.setStatus(restartableWorkflow.id, runningCollectors map { _.toKey }, ExecutionStatus.Starting)
          } yield ()
        }
      }
    }

    for {
      // Find all executions for the specified workflow that are not NotStarted or Done.
      executions <- globalDataAccess.getExecutionsForRestart(restartableWorkflow.id)
      // Examine statuses/types of executions, reset statuses as necessary.
      _ <- handleExecutionStatuses(executions)
    } yield ()
  }

  override def backendType = BackendType.JES

  def gcsAuthFilePath(descriptor: WorkflowDescriptor): String = {
    // If we are going to upload an auth file we need a valid GCS path passed via workflow options.
    val bucket = descriptor.workflowOptions.get(AuthFilePathOptionKey) getOrElse {
      throw new UnsatisfiedInputsException(s"$AuthFilePathOptionKey has not been found in workflow options.")
    }
    s"$bucket/${descriptor.id}_auth.json"
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
    descriptor.workflowOptions.getOrElse(GoogleProjectOptionKey, jesConf.project)
  }

  // Create an input parameter containing the path to this authentication file, if needed
  def gcsAuthParameter(descriptor: WorkflowDescriptor): Option[JesInput] = {
    if ((jesConf.localizeWithRefreshToken || jesConf.isDockerAuthenticated) && descriptor.workflowOptions.get(AuthFilePathOptionKey).isSuccess)
      Option(authGcsCredentialsPath(gcsAuthFilePath(descriptor)))
    else None
  }

  override def findResumableExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, JobKey]] = {
    globalDataAccess.findResumableJesExecutions(id)
  }

  /**
    * Will try to instantiate a user authenticated IOInterface if possible.
    * If not, will default back to the Cromwell authenticated IOInterface.
    */
  override def ioInterface(workflowOptions: WorkflowOptions): GoogleCloudStorage = {
    val userCredentials = for {
      secrets <- jesConf.googleSecrets
      token <- workflowOptions.get(RefreshTokenOptionKey).toOption
    } yield GoogleCredentialFactory.forRefreshToken(secrets, token)
    userCredentials map { GcsFactory(jesConf.applicationName, _) } getOrElse jesCromwellConnection.storage
  }
}
