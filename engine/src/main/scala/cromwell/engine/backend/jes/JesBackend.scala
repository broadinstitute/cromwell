package cromwell.engine.backend.jes

import java.net.SocketTimeoutException
import java.nio.file.{FileSystem, Path, Paths}

import akka.actor.ActorSystem
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.api.client.http.HttpResponseException
import com.google.api.client.util.ExponentialBackOff.Builder
import com.google.api.services.genomics.model.{LocalCopy, PipelineParameter}
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.{CallOutput, CallOutputs, WorkflowId, WorkflowOptions}
import cromwell.engine.ExecutionIndex.IndexEnhancedInt
import cromwell.engine.ExecutionStatus._
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.jes.Run.{RunStatus, TerminalRunStatus}
import cromwell.engine.backend.jes.authentication._
import cromwell.engine.db.DataAccess.{ExecutionKeyToJobKey, globalDataAccess}
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.db.slick.{Execution, ExecutionInfo}
import cromwell.engine.io.gcs._
import cromwell.filesystems.gcs.{GoogleConfiguration, RefreshTokenMode, GoogleAuthMode}
import cromwell.filesystems.gcs.RefreshTokenMode
import cromwell.logging.WorkflowLogger
import cromwell.util.{CromwellAggregatedException, SimpleExponentialBackoff, TryUtil}
import spray.json.JsObject
import wdl4s.AstTools.EnhancedAstNode
import wdl4s.command.ParameterCommandPart
import wdl4s.expression.NoFunctions
import wdl4s.values._
import wdl4s.{CallInputs, UnsatisfiedInputsException, _}

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object JesBackend {
  val ExecParamName = "exec"
  val MonitoringParamName = "monitoring"
  val ExtraConfigParamName = "__extra_config_gcs_path"
  val JesExecScript = "exec.sh"
  val JesMonitoringScript = "monitoring.sh"
  val JesMonitoringLogFile = "monitoring.log"

  // Workflow options keys
  val GcsRootOptionKey = "jes_gcs_root"
  val MonitoringScriptOptionKey = "monitoring_script"
  val GoogleProjectOptionKey = "google_project"
  val AuthFilePathOptionKey = "auth_bucket"
  val WriteToCacheOptionKey = "write_to_cache"
  val ReadFromCacheOptionKey = "read_from_cache"
  val OptionKeys = Set(
    GoogleConfiguration.RefreshTokenOptionKey, GcsRootOptionKey, MonitoringScriptOptionKey, GoogleProjectOptionKey,
    AuthFilePathOptionKey, WriteToCacheOptionKey, ReadFromCacheOptionKey
  ) ++ WorkflowDescriptor.OptionKeys

  def authGcsCredentialsPath(gcsPath: String): JesInput = JesLiteralInput(ExtraConfigParamName, gcsPath)

  // Only executions in Running state with a recorded operation ID are resumable.
  private val IsResumable: (Execution, Seq[ExecutionInfo]) => Boolean = (e: Execution, eis: Seq[ExecutionInfo]) => {
    e.status.toExecutionStatus == ExecutionStatus.Running &&
      eis.exists(ei => ei.key == JesBackend.InfoKeys.JesRunId && ei.value.isDefined)
  }

  /** In this context transient means that a job is Running but doesn't have a Jes Run ID yet.
    * This could happen if a call starts Running but Cromwell is stopped before the request to JES has been made
    * (or the ID persisted to the DB) .*/
  private val IsTransient: (Execution, Seq[ExecutionInfo]) => Boolean = (e: Execution, eis: Seq[ExecutionInfo]) => {
    e.status.toExecutionStatus == ExecutionStatus.Running &&
      eis.exists(ei => ei.key == JesBackend.InfoKeys.JesRunId && ei.value.isEmpty)
  }

  case class JesBackendJobKey(jesRunId: String) extends BackendJobKey

  private val BuildBackendJobKey: (Execution, Seq[ExecutionInfo]) => BackendJobKey = (e: Execution, eis: Seq[ExecutionInfo]) => {
    JesBackendJobKey(eis.find(_.key == JesBackend.InfoKeys.JesRunId).get.value.get)
  }

  /**
   * Takes a path in GCS and comes up with a local path which is unique for the given GCS path
   *
   * @param gcsPath The input path
   * @return A path which is unique per input path
   */
  def localFilePathFromCloudStoragePath(gcsPath: GcsPath): Path = {
    Paths.get(gcsPath.bucket).resolve(gcsPath.objectName)
  }

  /**
   * Takes a single WdlValue and maps google cloud storage (GCS) paths into an appropriate local file path.
   * If the input is not a WdlFile, or the WdlFile is not a GCS path, the mapping is a noop.
   *
   * @param wdlValue the value of the input
   * @return a new FQN to WdlValue pair, with WdlFile paths modified if appropriate.
   */
  def gcsPathToLocal(wdlValue: WdlValue): WdlValue = {
    wdlValue match {
      case wdlFile: WdlFile =>
        GcsPath.parse(wdlFile.value) match {
          case Success(gcsPath) => WdlFile(localFilePathFromCloudStoragePath(gcsPath).toString, wdlFile.isGlob)
          case Failure(e) => wdlValue
        }
      case wdlArray: WdlArray => wdlArray map gcsPathToLocal
      case wdlMap: WdlMap => wdlMap map { case (k, v) => gcsPathToLocal(k) -> gcsPathToLocal(v) }
      case _ => wdlValue
    }
  }

  def isFatalJesException(t: Throwable): Boolean = t match {
    case e: HttpResponseException if e.getStatusCode == 403 => true
    case _ => false
  }

  def isTransientJesException(t: Throwable): Boolean = t match {
      // Quota exceeded
    case e: HttpResponseException if e.getStatusCode == 429 => true
    case _ => false
  }

  protected def withRetry[T](f: Option[T] => T, logger: WorkflowLogger, failureMessage: String) = {
    TryUtil.retryBlock(
      fn = f,
      retryLimit = Option(10),
      backoff = SimpleExponentialBackoff(5 seconds, 10 seconds, 1.1D),
      logger = logger,
      failMessage = Option(failureMessage),
      isFatal = isFatalJesException,
      isTransient = isTransientJesException
    )
  }

  sealed trait JesParameter {
    def name: String
    def toGooglePipelineParameter: PipelineParameter
    def toGoogleRunParameter: String
  }

  sealed trait JesInput extends JesParameter

  final case class JesFileInput(name: String, gcs: String, local: Path, mount: JesAttachedDisk) extends JesInput {
    def toGooglePipelineParameter = {
      new PipelineParameter().setName(name).setLocalCopy(
        new LocalCopy().setDisk(mount.name).setPath(local.toString)
      )
    }
    val toGoogleRunParameter: String = gcs
    def containerPath: Path = mount.mountPoint.resolve(local)
  }

  final case class JesLiteralInput(name: String, value: String) extends JesInput {
    def toGooglePipelineParameter = new PipelineParameter().setName(name)
    val toGoogleRunParameter: String = value
  }

  final case class JesFileOutput(name: String, gcs: String, local: Path, mount: JesAttachedDisk) extends JesParameter {
    def toGooglePipelineParameter = {
      new PipelineParameter().setName(name).setLocalCopy(
        new LocalCopy().setDisk(mount.name).setPath(local.toString)
      )
    }
    val toGoogleRunParameter: String = gcs
  }

  implicit class EnhancedExecution(val execution: Execution) extends AnyVal {
    import cromwell.engine.ExecutionIndex._
    def toKey: ExecutionDatabaseKey = ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex, execution.attempt)
    def isScatter: Boolean = execution.callFqn.contains(Scatter.FQNIdentifier)
    def executionStatus: ExecutionStatus = ExecutionStatus.withName(execution.status)
  }

  object InfoKeys {
    val JesRunId = "JES_RUN_ID"
    val JesStatus = "JES_STATUS"
  }

  def jesLogBasename(key: JobKey) = {
    val index = key.index.map(s => s"-$s").getOrElse("")
    s"${key.scope.unqualifiedName}$index"
  }

  def jesLogStdoutFilename(key: JobKey) = s"${jesLogBasename(key)}-stdout.log"
  def jesLogStderrFilename(key: JobKey) = s"${jesLogBasename(key)}-stderr.log"
  def jesLogFilename(key: JobKey) = s"${jesLogBasename(key)}.log"
  def jesReturnCodeFilename(key: JobKey) = s"${jesLogBasename(key)}-rc.txt"
  def globDirectory(glob: String) = s"glob-${glob.md5Sum}/"
  def rootPath(options: WorkflowOptions): String = options.getOrElse(GcsRootOptionKey, ProductionJesConfiguration.jesConf.executionBucket).stripSuffix("/")
  def globOutputPath(callPath: Path, glob: String) = callPath.resolve(s"glob-${glob.md5Sum}/")
  def preemptible(jobDescriptor: BackendCallJobDescriptor, maxPreemption: Int) = jobDescriptor.key.attempt <= maxPreemption
  def monitoringIO(jobDescriptor: BackendCallJobDescriptor): Option[JesInput] = {
    jobDescriptor.workflowDescriptor.workflowOptions.get(MonitoringScriptOptionKey) map { path =>
      JesFileInput(s"$MonitoringParamName-in", GcsPath(path).toString, Paths.get(JesMonitoringScript), jobDescriptor.workingDisk)
    } toOption
  }

  /**
    * Generates a json containing auth information based on the parameters provided.
 *
    * @return a string representation of the json
    */
  def generateAuthJson(authInformation: Option[JesAuthInformation]*) = {
    authInformation.flatten map { _.toMap } match {
      case Nil => None
      case jsons =>
        val authsValues = jsons.reduce(_ ++ _) mapValues JsObject.apply
        Option(JsObject("auths" -> JsObject(authsValues)).prettyPrint)
    }
  }
}

/**
 * Representing a running JES execution, instances of this class are never Done and it is never okay to
 * ask them for results.
 */
case class JesPendingExecutionHandle(jobDescriptor: BackendCallJobDescriptor,
                                     jesOutputs: Seq[JesFileOutput],
                                     run: Run,
                                     previousStatus: Option[RunStatus]) extends ExecutionHandle {
  override val isDone = false
  override val result = NonRetryableExecution(new IllegalStateException("JesPendingExecutionHandle cannot yield a result"))
}

case class JesBackend(actorSystem: ActorSystem)
  extends Backend
  with LazyLogging
  with ProductionJesAuthentication
  with ProductionJesConfiguration {

  import backend.io._
  import better.files._

  /**
    * Exponential Backoff Builder to be used when polling for call status.
    */
  final private lazy val pollBackoffBuilder = new Builder()
    .setInitialIntervalMillis(20.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Int.MaxValue)
    .setMaxIntervalMillis(10.minutes.toMillis.toInt)
    .setMultiplier(1.1D)

  override def pollBackoff = pollBackoffBuilder.build()

  override def rootPath(options: WorkflowOptions) = options.getOrElse(GcsRootOptionKey, ProductionJesConfiguration.jesConf.executionBucket).stripSuffix("/")

  // FIXME: Add proper validation of jesConf and have it happen up front to provide fail-fast behavior (will do as a separate PR)

  override def adjustInputPaths(jobDescriptor: BackendCallJobDescriptor): CallInputs = jobDescriptor.locallyQualifiedInputs mapValues gcsPathToLocal

  private def writeAuthenticationFile(workflow: WorkflowDescriptor): Try[Unit] = {
    val log = workflowLogger(workflow)

    generateAuthJson(dockerConf, getGcsAuthInformation(workflow)) map { content =>

      val path = gcsAuthFilePath(workflow)
      def upload(prev: Option[Unit]): Unit = path.writeAsJson(content)

      log.info(s"Creating authentication file for workflow ${workflow.id} at \n ${path.toString}")
      withRetry(upload, log, s"Exception occurred while uploading auth file to $path")
    } getOrElse Success(())
  }

  def engineFunctions(fileSystems: List[FileSystem], workflowContext: WorkflowContext): WorkflowEngineFunctions = {
    new JesWorkflowEngineFunctions(fileSystems, workflowContext)
  }

  /**
   * Get a GcsLocalizing from workflow options if client secrets and refresh token are available.
   */
  def getGcsAuthInformation(workflow: WorkflowDescriptor): Option[JesAuthInformation] = {
    def extractSecrets(gcsMode: GoogleAuthMode) = gcsMode match {
      case RefreshTokenMode(secrets) => Option(secrets)
      case _ => None
    }

    for {
      secrets <- extractSecrets(gcsConf.authMode)
      token <- workflow.workflowOptions.get(GoogleConfiguration.RefreshTokenOptionKey).toOption
    } yield GcsLocalizing(secrets, token)
  }

  /*
   * No need to copy GCS inputs for the workflow we should be able to directly reference them
   * Create an authentication json file containing docker credentials and/or user account information
   */
  override def initializeForWorkflow(workflow: WorkflowDescriptor): Try[Unit] = {
    writeAuthenticationFile(workflow)
  }

  override def assertWorkflowOptions(options: WorkflowOptions): Unit = {
    // Warn for unrecognized option keys
    options.toMap.keySet.diff(OptionKeys) match {
      case unknowns if unknowns.nonEmpty => logger.warn(s"Unrecognized workflow option(s): ${unknowns.mkString(", ")}")
      case _ =>
    }

    gcsConf.authMode match {
      case _: RefreshTokenMode => Seq(GoogleConfiguration.RefreshTokenOptionKey) filterNot options.toMap.keySet match {
        case missing if missing.nonEmpty =>
          throw new IllegalArgumentException(s"Missing parameters in workflow options: ${missing.mkString(", ")}")
        case _ =>
      }
      case _ =>
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

  private def deleteAuthFile(authFilePath: Path, log: WorkflowLogger): Future[Unit] = {
      def gcsCheckAuthFileExists(prior: Option[Boolean]): Boolean = authFilePath.exists
      def gcsAttemptToDeleteObject(prior: Option[Unit]): Unit = authFilePath.delete()
      withRetry(gcsCheckAuthFileExists, log, s"Failed to query for auth file: $authFilePath") match {
        case Success(exists) if exists =>
          withRetry(gcsAttemptToDeleteObject, log, s"Failed to delete auth file: $authFilePath") match {
            case Success(_) => Future.successful(Unit)
            case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 =>
              log.warn(s"Could not delete the auth file $authFilePath: File Not Found")
              Future.successful(Unit)
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

  def stdoutStderr(jobDescriptor: BackendCallJobDescriptor): CallLogs = {
    CallLogs(
      stdout = WdlFile(jobDescriptor.jesStdoutGcsPath.toString),
      stderr = WdlFile(jobDescriptor.jesStderrGcsPath.toString),
      Option(Map("log" -> WdlFile(jobDescriptor.jesLogGcsPath.toString)))
    )
  }

  private def executeOrResume(jobDescriptor: BackendCallJobDescriptor, runIdForResumption: Option[String])(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    val log = jobLogger(jobDescriptor)
    log.info(s"Call GCS path: ${jobDescriptor.callRootPath}")
    val monitoringScript: Option[JesInput] = monitoringIO(jobDescriptor)
    val monitoringOutput = monitoringScript map { _ =>
      JesFileOutput(s"$MonitoringParamName-out", jobDescriptor.defaultMonitoringOutputPath.toString, Paths.get(JesMonitoringLogFile), jobDescriptor.workingDisk)
    }

    val jesInputs: Seq[JesInput] = generateJesInputs(jobDescriptor).toSeq ++ monitoringScript :+ jobDescriptor.cmdInput
    val jesOutputs: Seq[JesFileOutput] = generateJesOutputs(jobDescriptor) ++ monitoringOutput

    jobDescriptor.instantiateCommand match {
      case Success(command) => runWithJes(jobDescriptor, command, jesInputs, jesOutputs, runIdForResumption, monitoringScript.isDefined)
      case Failure(ex: SocketTimeoutException) => throw ex // probably a GCS transient issue, throwing will cause it to be retried
      case Failure(ex) => FailedExecutionHandle(ex)
    }
  }

  def execute(jobDescriptor: BackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[ExecutionHandle] = executeOrResume(jobDescriptor, runIdForResumption = None)

  def resume(jobDescriptor: BackendCallJobDescriptor, jobKey: BackendJobKey)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    val runId = Option(jobKey) collect { case jesKey: JesBackendJobKey => jesKey.jesRunId }
    executeOrResume(jobDescriptor, runIdForResumption = runId)
  }

  def useCachedCall(cachedCallDescriptor: BackendCallJobDescriptor, callDescriptor: BackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {

    import better.files._
    // FIXME Temporary workaround so both paths share the same FileSystemProvider, which allows copy in the cloud instead of streaming through cromwell
    val jobDescriptorPath = callRootPath(callDescriptor).toString.toDirectory(cachedCallDescriptor.workflowDescriptor.fileSystems)

    def renameCallSpecificFiles = {
      val namesBuilder = List(
        JesBackend.jesLogStdoutFilename _,
        JesBackend.jesLogStderrFilename _,
        JesBackend.jesLogFilename _,
        JesBackend.jesReturnCodeFilename _
      )

      /* Paths of the files in the cached call. */
      val cachedPaths = namesBuilder map { _(cachedCallDescriptor.key) } map callRootPath(cachedCallDescriptor).resolve
      /* Expected names of the files in the current call.
       * They might be different if the cached call doesn't have the same name / shard
       * as the current call, because this information is part of the filename (TODO should it ?).
       */
      val backendNames = namesBuilder map { _(callDescriptor.key) }

      (cachedPaths zip backendNames) map {
        case (cachedPath, backendName) => cachedPath.renameTo(backendName)
      }
    }

    def copyingWork: Try[Unit] = for {
        _ <- Try(callRootPath(cachedCallDescriptor).copyTo(jobDescriptorPath))
        _ <- Try(renameCallSpecificFiles)
    } yield ()

    Future {
      val log = jobLogger(callDescriptor)
        copyingWork match {
          case Failure(ex) =>
            log.error(s"Exception occurred while attempting to copy outputs from ${cachedCallDescriptor.callRootPath} to $jobDescriptorPath", ex)
            FailedExecutionHandle(ex).future
          case Success(_) => postProcess(callDescriptor) match {
            case Success(outputs) => callDescriptor.hash map { h =>
              SuccessfulExecutionHandle(outputs, List.empty[ExecutionEventEntry], cachedCallDescriptor.downloadRcFile.get.stripLineEnd.toInt, h, Option(cachedCallDescriptor)) }
            case Failure(ex: CromwellAggregatedException) if ex.throwables collectFirst { case s: SocketTimeoutException => s } isDefined =>
              // TODO: What can we return here to retry this operation?
              // TODO: This match clause is similar to handleSuccess(), though it's subtly different for this specific case
              val error = "Socket timeout occurred in evaluating one or more of the output expressions"
              log.error(error, ex)
              FailedExecutionHandle(new Exception(error, ex)).future
            case Failure(ex) => FailedExecutionHandle(ex).future
          }
      }
    } flatten
  }

  /**
    * Turns WdlFiles into relative paths.  These paths are relative to the working disk
    *
    * relativeLocalizationPath("foo/bar.txt") -> "foo/bar.txt"
    * relativeLocalizationPath("gs://some/bucket/foo.txt") -> "some/bucket/foo.txt"
    */
  private def relativeLocalizationPath(file: WdlFile): WdlFile = {
    GcsPath.parse(file.value) match {
      case Success(gcsPath) => WdlFile(gcsPath.bucket + "/" + gcsPath.objectName, file.isGlob)
      case Failure(e) => file
    }
  }

  def generateJesInputs(jobDescriptor: BackendCallJobDescriptor): Iterable[JesInput] = {
    /**
      * Commands in WDL tasks can also generate input files.  For example: ./my_exec --file=${write_lines(arr)}
      *
      * write_lines(arr) would produce a string-ified version of the array stored as a GCS path.  The next block of code
      * will go through each ${...} expression within the task's command section and find all write_*() ASTs and
      * evaluate them so the files are written to GCS and the they can be included as inputs to Google's Pipeline object
      */
    val commandExpressions = jobDescriptor.key.scope.task.commandTemplate.collect({
      case x: ParameterCommandPart => x.expression
    })

    val writeFunctionAsts = commandExpressions.map(_.ast).flatMap(x => AstTools.findAsts(x, "FunctionCall")).collect({
      case y if y.getAttribute("name").sourceString.startsWith("write_") => y
    })

    val evaluatedExpressionMap = writeFunctionAsts map { ast =>
      val expression = WdlExpression(ast)
      val value = expression.evaluate(jobDescriptor.lookupFunction(Map.empty), jobDescriptor.callEngineFunctions)
      expression.toWdlString.md5SumShort -> value
    } toMap

    val writeFunctionFiles = evaluatedExpressionMap collect { case (k, v: Success[_]) => k -> v.get } collect { case (k, v: WdlFile) => k -> Seq(v)}

    /** Collect all WdlFiles from inputs to the call */
    val callInputFiles = jobDescriptor.locallyQualifiedInputs mapValues { _.collectAsSeq { case w: WdlFile => w } }

    (callInputFiles ++ writeFunctionFiles) flatMap {
      case (name, files) => jesInputsFromWdlFiles(name, files, files.map(relativeLocalizationPath), jobDescriptor)
    }
  }

  /**
   * Takes two arrays of remote and local WDL File paths and generates the necessary JesInputs.
   */
  private def jesInputsFromWdlFiles(jesNamePrefix: String,
                                    remotePathArray: Seq[WdlFile],
                                    localPathArray: Seq[WdlFile],
                                    jobDescriptor: BackendCallJobDescriptor): Iterable[JesInput] = {
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((remotePath, localPath), index) =>
        Seq(JesFileInput(s"$jesNamePrefix-$index", remotePath.valueString, Paths.get(localPath.valueString), jobDescriptor.workingDisk))
    }
  }

  def generateJesOutputs(jobDescriptor: BackendCallJobDescriptor): Seq[JesFileOutput] = {
    val log = jobLogger(jobDescriptor)
    val wdlFileOutputs = jobDescriptor.key.scope.task.outputs flatMap { taskOutput =>
      taskOutput.requiredExpression.evaluateFiles(jobDescriptor.lookupFunction(Map.empty), NoFunctions, taskOutput.wdlType) match {
        case Success(wdlFiles) => wdlFiles map relativeLocalizationPath
        case Failure(ex) =>
          log.warn(s"Could not evaluate $taskOutput: ${ex.getMessage}", ex)
          Seq.empty[WdlFile]
      }
    }

    // Create the mappings. GLOB mappings require special treatment (i.e. stick everything matching the glob in a folder)
    wdlFileOutputs.distinct map { wdlFile =>
      val destination = wdlFile match {
        case WdlSingleFile(filePath) => jobDescriptor.callRootPath.resolve(filePath).toString
        case WdlGlobFile(filePath) => jobDescriptor.globOutputPath(filePath)
      }
      val (relpath, disk) = relativePathAndAttachedDisk(wdlFile.value, jobDescriptor.callRuntimeAttributes.disks)
      JesFileOutput(makeSafeJesReferenceName(wdlFile.value), destination, relpath, disk)
    }
  }

  /**
    * Given a path (relative or absolute), returns a (Path, JesAttachedDisk) tuple where the Path is
    * relative to the AttachedDisk's mount point
    *
    * @throws Exception if the `path` does not live in one of the supplied `disks`
    */
  private def relativePathAndAttachedDisk(path: String, disks: Seq[JesAttachedDisk]): (Path, JesAttachedDisk) = {
    val absolutePath = Paths.get(path) match {
      case p if !p.isAbsolute => Paths.get(JesWorkingDisk.MountPoint).resolve(p)
      case p => p
    }

    disks.find(d => absolutePath.startsWith(d.mountPoint)) match {
      case Some(disk) => (disk.mountPoint.relativize(absolutePath), disk)
      case None =>
        throw new Exception(s"Absolute path $path doesn't appear to be under any mount points: ${disks.map(_.toString).mkString(", ")}")
    }
  }

  /**
   * If the desired reference name is too long, we don't want to break JES or risk collisions by arbitrary truncation. So,
   * just use a hash. We only do this when needed to give better traceability in the normal case.
   */
  private def makeSafeJesReferenceName(referenceName: String) = {
    if (referenceName.length <= 127) referenceName else referenceName.md5Sum
  }

  private def uploadCommandScript(jobDescriptor: BackendCallJobDescriptor, command: String, withMonitoring: Boolean): Try[Unit] = {
    val monitoring = if (withMonitoring) {
      s"""|touch $JesMonitoringLogFile
          |chmod u+x $JesMonitoringScript
          |./$JesMonitoringScript > $JesMonitoringLogFile &""".stripMargin
    } else ""

    val tmpDir = Paths.get(JesWorkingDisk.MountPoint).resolve("tmp")
    val rcPath = Paths.get(JesWorkingDisk.MountPoint).resolve(JesBackend.jesReturnCodeFilename(jobDescriptor.key))

    val fileContent =
      s"""
         |#!/bin/bash
         |export _JAVA_OPTIONS=-Djava.io.tmpdir=$tmpDir
         |export TMPDIR=$tmpDir
         |cd ${JesWorkingDisk.MountPoint}
         |$monitoring
         |$command
         |echo $$? > $rcPath
       """.stripMargin.trim

    def attemptToUploadObject(priorAttempt: Option[Unit]): Unit = Option(jobDescriptor.gcsExecPath.write(fileContent))

    val log = jobLogger(jobDescriptor)
    withRetry(attemptToUploadObject, log, s"${jobLogger(jobDescriptor).tag} Exception occurred while uploading script to ${jobDescriptor.gcsExecPath}")
  }

  private def createJesRun(jobDescriptor: BackendCallJobDescriptor, jesParameters: Seq[JesParameter], runIdForResumption: Option[String]): Try[Run] = {
      def attemptToCreateJesRun(priorAttempt: Option[Run]): Run = Pipeline(
        jobDescriptor,
        jesParameters,
        googleProject(jobDescriptor.workflowDescriptor),
        genomicsInterface,
        runIdForResumption
      ).run

      val log = jobLogger(jobDescriptor)
      if (jobDescriptor.preemptible) log.info("Starting call with pre-emptible VM")
      withRetry(attemptToCreateJesRun, log, "Exception occurred while creating JES Run")
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
  private def gcsInputToGcsOutput(jobDescriptor: BackendCallJobDescriptor, inputValue: WdlValue): WdlValue = {
    // Convert to the local path where the file is localized to in the VM:
    val vmLocalizationPath = gcsPathToLocal(inputValue)

    vmLocalizationPath match {
      // If it's a file, work out where the file would be delocalized to, otherwise no-op:
      case x : WdlFile =>
        val delocalizationPath = callRootPath(jobDescriptor).resolve(vmLocalizationPath.valueString).toString
        WdlFile(delocalizationPath)
      case a: WdlArray => WdlArray(a.wdlType, a.value map { f => gcsInputToGcsOutput(jobDescriptor, f) })
      case m: WdlMap => WdlMap(m.wdlType, m.value map { case (k, v) => gcsInputToGcsOutput(jobDescriptor, k) -> gcsInputToGcsOutput(jobDescriptor, v) })
      case other => other
    }
  }

  private def customLookupFunction(jobDescriptor: BackendCallJobDescriptor, alreadyGeneratedOutputs: Map[String, WdlValue]): String => WdlValue = toBeLookedUp => {
    val originalLookup = jobDescriptor.lookupFunction(alreadyGeneratedOutputs)
    gcsInputToGcsOutput(jobDescriptor, originalLookup(toBeLookedUp))
  }

  def wdlValueToGcsPath(jesOutputs: Seq[JesFileOutput])(value: WdlValue): WdlValue = {
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

  def postProcess(jobDescriptor: BackendCallJobDescriptor): Try[CallOutputs] = {
    val outputs = jobDescriptor.call.task.outputs
    val outputFoldingFunction = getOutputFoldingFunction(jobDescriptor)
    val outputMappings = outputs.foldLeft(Seq.empty[AttemptedLookupResult])(outputFoldingFunction).map(_.toPair).toMap
    TryUtil.sequenceMap(outputMappings) map { outputMap =>
      outputMap mapValues { v =>
        CallOutput(v, jobDescriptor.workflowDescriptor.hash(v))
      }
    }
  }

  private def getOutputFoldingFunction(jobDescriptor: BackendCallJobDescriptor): (Seq[AttemptedLookupResult], TaskOutput) => Seq[AttemptedLookupResult] = {
    (currentList: Seq[AttemptedLookupResult], taskOutput: TaskOutput) => {
      currentList ++ Seq(AttemptedLookupResult(taskOutput.name, outputLookup(taskOutput, jobDescriptor, currentList)))
    }
  }

  private def outputLookup(taskOutput: TaskOutput, jobDescriptor: BackendCallJobDescriptor, currentList: Seq[AttemptedLookupResult]) = for {
  /**
    * This will evaluate the task output expression and coerces it to the task output's type.
    * If the result is a WdlFile, then attempt to find the JesOutput with the same path and
    * return a WdlFile that represents the GCS path and not the local path.  For example,
    *
    * <pre>
    * output {
    *   File x = "out" + ".txt"
    * }
    * </pre>
    *
    * "out" + ".txt" is evaluated to WdlString("out.txt") and then coerced into a WdlFile("out.txt")
    * Then, via wdlFileToGcsPath(), we attempt to find the JesOutput with .name == "out.txt".
    * If it is found, then WdlFile("gs://some_bucket/out.txt") will be returned.
    */
    wdlValue <- taskOutput.requiredExpression.evaluate(customLookupFunction(jobDescriptor, currentList.toLookupMap), jobDescriptor.callEngineFunctions)
    coercedValue <- taskOutput.wdlType.coerceRawValue(wdlValue)
    value = wdlValueToGcsPath(generateJesOutputs(jobDescriptor))(coercedValue)
  } yield value

  def executionResult(status: RunStatus, handle: JesPendingExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    val log = jobLogger(handle.jobDescriptor)

    try {
      val jobDescriptor = handle.jobDescriptor
      val outputMappings = postProcess(jobDescriptor)
      lazy val stderrLength: Long = jobDescriptor.jesStderrGcsPath.size
      lazy val returnCodeContents = jobDescriptor.downloadRcFile
      lazy val returnCode = returnCodeContents map { _.trim.toInt }
      lazy val continueOnReturnCode = jobDescriptor.callRuntimeAttributes.continueOnReturnCode

      status match {
        case Run.Success(events) if jobDescriptor.callRuntimeAttributes.failOnStderr && stderrLength.intValue > 0 =>
          // returnCode will be None if it couldn't be downloaded/parsed, which will yield a null in the DB
          FailedExecutionHandle(new Throwable(s"${log.tag} execution failed: stderr has length $stderrLength"), returnCode.toOption).future
        case Run.Success(events) if returnCodeContents.isFailure =>
          val exception = returnCode.failed.get
          log.warn(s"${log.tag} could not download return code file, retrying: " + exception.getMessage, exception)
          // Return handle to try again.
          handle.future
        case Run.Success(events) if returnCode.isFailure =>
          FailedExecutionHandle(new Throwable(s"${log.tag} execution failed: could not parse return code as integer: " + returnCodeContents.get)).future
        case Run.Success(events) if !continueOnReturnCode.continueFor(returnCode.get) =>
          FailedExecutionHandle(new Throwable(s"${log.tag} execution failed: disallowed command return code: " + returnCode.get), returnCode.toOption).future
        case Run.Success(events) =>
          jobDescriptor.hash map { h => handleSuccess(outputMappings, jobDescriptor.workflowDescriptor, events, returnCode.get, h, handle) }
        case Run.Failed(errorCode, errorMessage, events) => handleFailure(jobDescriptor, errorCode, errorMessage, events, log)
      }
    } catch {
      case e: Exception =>
        log.warn("Caught exception trying to download result, retrying: " + e.getMessage, e)
        // Return the original handle to try again.
        handle.future
    }
  } flatten

  private def runWithJes(jobDescriptor: BackendCallJobDescriptor,
                         command: String,
                         jesInputs: Seq[JesInput],
                         jesOutputs: Seq[JesFileOutput],
                         runIdForResumption: Option[String],
                         withMonitoring: Boolean): ExecutionHandle = {
    val log = jobLogger(jobDescriptor)
    val jesParameters = jobDescriptor.standardParameters ++ gcsAuthParameter(jobDescriptor.workflowDescriptor) ++ jesInputs ++ jesOutputs
    log.info(s"`$command`")

    val jesJobSetup = for {
      _ <- uploadCommandScript(jobDescriptor, command, withMonitoring)
      run <- createJesRun(jobDescriptor, jesParameters, runIdForResumption)
    } yield run

    jesJobSetup match {
      case Failure(ex) =>
        log.warn(s"Failed to create a JES run", ex)
        throw ex  // Probably a transient issue, throwing retries it
      case Success(run) => JesPendingExecutionHandle(jobDescriptor, jesOutputs, run, previousStatus = None)
    }
  }

  private def handleSuccess(outputMappings: Try[CallOutputs],
                            workflowDescriptor: WorkflowDescriptor,
                            executionEvents: Seq[ExecutionEventEntry],
                            returnCode: Int,
                            hash: ExecutionHash,
                            executionHandle: ExecutionHandle): ExecutionHandle = {
    outputMappings match {
      case Success(outputs) => SuccessfulExecutionHandle(outputs, executionEvents, returnCode, hash)
      case Failure(ex: CromwellAggregatedException) if ex.throwables collectFirst { case s: SocketTimeoutException => s } isDefined =>
        // Return the execution handle in this case to retry the operation
        executionHandle
      case Failure(ex) => FailedExecutionHandle(ex)
    }
  }

  private def extractErrorCodeFromErrorMessage(errorMessage: String): Int = {
    errorMessage.substring(0, errorMessage.indexOf(':')).toInt
  }

  private def preempted(errorCode: Int, errorMessage: Option[String], jobDescriptor: BackendCallJobDescriptor, logger: WorkflowLogger): Boolean = {
    def isPreemptionCode(code: Int) = code == 13 || code == 14

    try {
      errorCode == 10 && errorMessage.isDefined && isPreemptionCode(extractErrorCodeFromErrorMessage(errorMessage.get)) && jobDescriptor.preemptible
    } catch {
      case _: NumberFormatException | _: StringIndexOutOfBoundsException =>
        logger.warn(s"Unable to parse JES error code from error message: ${errorMessage.get}, assuming this was not a preempted VM.")
        false
    }
  }

  private def handleFailure(jobDescriptor: BackendCallJobDescriptor, errorCode: Int, errorMessage: Option[String], events: Seq[ExecutionEventEntry], logger: WorkflowLogger) = {
    import lenthall.numeric.IntegerUtil._

    val taskName = s"${jobDescriptor.workflowDescriptor.id}:${jobDescriptor.call.unqualifiedName}"
    val attempt = jobDescriptor.key.attempt

    if (errorMessage.exists(_.contains("Operation canceled at")))  {
      AbortedExecutionHandle.future
    } else if (preempted(errorCode, errorMessage, jobDescriptor, logger)) {
      val preemptedMsg = s"Task $taskName was preempted for the ${attempt.toOrdinal} time."
      val max = jobDescriptor.maxPreemption

      if (attempt < max) {
        val e = new PreemptedException(
          s"""$preemptedMsg The call will be restarted with another preemptible VM (max preemptible attempts number is $max).
             |Error code $errorCode. Message: $errorMessage""".stripMargin
        )
        RetryableExecutionHandle(e, None, events).future
      } else {
        val e = new PreemptedException(
          s"""$preemptedMsg The maximum number of preemptible attempts ($max) has been reached. The call will be restarted with a non-preemptible VM.
             |Error code $errorCode. Message: $errorMessage)""".stripMargin)
        RetryableExecutionHandle(e, None, events).future
      }
    } else {
      val e = new Throwable(s"Task ${jobDescriptor.workflowDescriptor.id}:${jobDescriptor.call.unqualifiedName} failed: error code $errorCode. Message: ${errorMessage.getOrElse("null")}")
      FailedExecutionHandle(e, None, events).future
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
    import cromwell.engine.backend.jes.JesBackend.EnhancedExecution

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
            _ <- globalDataAccess.resetTransientExecutions(restartableWorkflow.id, IsTransient)
            _ <- globalDataAccess.setStartingStatus(restartableWorkflow.id, runningCollectors map { _.toKey })
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

  def gcsAuthFilePath(descriptor: WorkflowDescriptor): Path =  {
    // If we are going to upload an auth file we need a valid GCS path passed via workflow options.
    val bucket = descriptor.workflowOptions.get(AuthFilePathOptionKey) getOrElse descriptor.workflowRootPath.toString
    bucket.toPath(cromwellGcsFileSystem).resolve(s"${descriptor.id}_auth.json")
  }

  def googleProject(descriptor: WorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse(GoogleProjectOptionKey, jesConf.project)
  }

  // Create an input parameter containing the path to this authentication file, if needed
  def gcsAuthParameter(descriptor: WorkflowDescriptor): Option[JesInput] = {
    if (gcsConf.authMode.isInstanceOf[RefreshTokenMode] || dockerConf.isDefined)
      Option(authGcsCredentialsPath(gcsAuthFilePath(descriptor).toString))
    else None
  }

  override def findResumableExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[ExecutionKeyToJobKey]] = {
    globalDataAccess.findResumableExecutions(id, IsResumable, BuildBackendJobKey)
  }

  override def executionInfoKeys: List[String] = List(JesBackend.InfoKeys.JesRunId, JesBackend.InfoKeys.JesStatus)

  override def callEngineFunctions(descriptor: BackendCallJobDescriptor): CallEngineFunctions = {
    lazy val callRootPath = descriptor.callRootPath.toString
    lazy val jesStdoutPath = descriptor.jesStdoutGcsPath.toString
    lazy val jesStderrPath = descriptor.jesStderrGcsPath.toString

    val callContext = new CallContext(callRootPath, jesStdoutPath, jesStderrPath)
    new JesCallEngineFunctions(descriptor.workflowDescriptor.fileSystems, callContext)
  }

  override def fileSystems(options: WorkflowOptions): List[FileSystem] = List(userGcsFileSystem(options))

  override def instantiateCommand(descriptor: BackendCallJobDescriptor): Try[String] = {
    val backendInputs = adjustInputPaths(descriptor)
    descriptor.call.instantiateCommandLine(backendInputs, descriptor.callEngineFunctions, JesBackend.gcsPathToLocal)
  }

  override def poll(jobDescriptor: BackendCallJobDescriptor, previous: ExecutionHandle)(implicit ec: ExecutionContext) = Future {
    previous match {
      case handle: JesPendingExecutionHandle =>
        val wfId = handle.jobDescriptor.workflowDescriptor.shortId
        val tag = handle.jobDescriptor.key.tag
        val runId = handle.run.runId
        logger.debug(s"[UUID($wfId)$tag] Polling JES Job $runId")
        val status = Try(handle.run.checkStatus(jobDescriptor, handle.previousStatus))
        status match {
          case Success(s: TerminalRunStatus) => executionResult(s, handle)
          case Success(s) => handle.copy(previousStatus = Option(s)).future // Copy the current handle with updated previous status.
          case Failure(e: GoogleJsonResponseException) if e.getStatusCode == 404 =>
            logger.error(s"JES Job ID ${handle.run.runId} has not been found, failing call")
            FailedExecutionHandle(e).future
          case Failure(e: Exception) =>
            // Log exceptions and return the original handle to try again.
            logger.warn("Caught exception, retrying: " + e.getMessage, e)
            handle.future
          case Failure(e: Error) => Future.failed(e) // JVM-ending calamity.
          case Failure(throwable) =>
            // Someone has subclassed Throwable directly?
            FailedExecutionHandle(throwable).future
        }
      case f: FailedExecutionHandle => f.future
      case s: SuccessfulExecutionHandle => s.future
      case badHandle => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $badHandle"))
    }
  } flatten
}
