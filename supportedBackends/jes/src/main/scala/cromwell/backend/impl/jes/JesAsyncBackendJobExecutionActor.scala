package cromwell.backend.impl.jes


import java.net.SocketTimeoutException
import java.nio.file.{Path, Paths}
import java.util.UUID

import akka.actor.{Actor, ActorLogging, Props}
import better.files._
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.api.client.http.HttpResponseException
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.async.AsyncBackendJobExecutionActor.ExecutionMode
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, FailedRetryableExecutionHandle, NonRetryableExecution, SuccessfulExecutionHandle}
import cromwell.backend.impl.jes.RunStatus
import cromwell.backend.impl.jes.RunStatus.TerminalRunStatus

import cromwell.backend.impl.jes.authentication.JesDockerCredentials
import cromwell.backend.impl.jes.JesImplicits.{GoogleAuthWorkflowOptions, PathString}
import cromwell.backend.impl.jes.io._
import cromwell.backend.{AttemptedLookupResult, BackendConfigurationDescriptor, BackendJobDescriptor, BackendWorkflowDescriptor, ExecutionHash, PreemptedException}
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}
import cromwell.core.{CallOutput, CromwellAggregatedException, _}
import cromwell.filesystems.gcs.GoogleConfiguration
import wdl4s.AstTools._
import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s._
import wdl4s.command.ParameterCommandPart
import wdl4s.expression.NoFunctions
import wdl4s.util.TryUtil
import wdl4s.values._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object JesAsyncBackendJobExecutionActor {

  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor, completionPromise: Promise[BackendJobExecutionResponse]): Props = {
    Props(new JesAsyncBackendJobExecutionActor(jobDescriptor, configurationDescriptor, completionPromise))
  }

  object WorkflowOptionKeys {
    val MonitoringScript = "monitoring_script"
    val AuthFilePath = "auth_bucket"
    val GoogleProject = "google_project"
  }

  private val ExecParamName = "exec"
  private val MonitoringParamName = "monitoring"

  private val JesMonitoringScript = "monitoring.sh"
  private val JesMonitoringLogFile = "monitoring.log"
  private val JesExecScript = "exec.sh"
  private val ExtraConfigParamName = "__extra_config_gcs_path"

  private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
  }

  /**
    * Takes a path in GCS and comes up with a local path which is unique for the given GCS path
    *
    * @param gcsPath The input path
    * @return A path which is unique per input path
    */
  private def localFilePathFromCloudStoragePath(gcsPath: GcsPath): Path = {
    Paths.get(gcsPath.bucket).resolve(gcsPath.objectName)
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
        GcsPath.parse(wdlFile.value) match {
          case Success(gcsPath) => WdlFile(localFilePathFromCloudStoragePath(gcsPath).toString, wdlFile.isGlob)
          case Failure(e) => wdlValue
        }
      case wdlArray: WdlArray => wdlArray map gcsPathToLocal
      case wdlMap: WdlMap => wdlMap map { case (k, v) => gcsPathToLocal(k) -> gcsPathToLocal(v) }
      case _ => wdlValue
    }
  }

  /**
    * Representing a running JES execution, instances of this class are never Done and it is never okay to
    * ask them for results.
    */
  case class JesPendingExecutionHandle(jobDescriptor: BackendJobDescriptor,
                                       jesOutputs: Seq[JesFileOutput],
                                       run: Run,
                                       previousStatus: Option[RunStatus]) extends ExecutionHandle {
    override val isDone = false
    override val result = NonRetryableExecution(new IllegalStateException("JesPendingExecutionHandle cannot yield a result"))
  }
}


class JesAsyncBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                       val configurationDescriptor: BackendConfigurationDescriptor,
                                       override val completionPromise: Promise[BackendJobExecutionResponse])
  extends Actor with ActorLogging with AsyncBackendJobExecutionActor {

  import JesAsyncBackendJobExecutionActor._
  /**
    * Exponential Backoff Builder to be used when polling for call status.
    */
  override lazy val backoff = SimpleExponentialBackoff(initialInterval = 30 seconds, maxInterval = 60 seconds, multiplier = 1.1)

  private lazy val monitoringScript: Option[JesInput] = {
    jobDescriptor.descriptor.workflowOptions.get(WorkflowOptionKeys.MonitoringScript) map { path =>
      JesFileInput(s"$MonitoringParamName-in", GcsPath(path).toString, Paths.get(JesMonitoringScript), workingDisk)
    } toOption
  }
  private lazy val workflowDescriptor = jobDescriptor.descriptor
  private lazy val gcsFileSystem = buildGcsFileSystem(configurationDescriptor, jobDescriptor.descriptor)
  private lazy val backendConfig = configurationDescriptor.backendConfig
  private lazy val jesCallPaths = JesCallPaths(jobDescriptor.key, gcsFileSystem, jobDescriptor.descriptor, configurationDescriptor.backendConfig)
  private lazy val googleConfig = GoogleConfiguration(configurationDescriptor.globalConfig)
  private lazy val genomicsFactory = GenomicsFactory(googleConfig, jesAttributes.genomicsAuth.credential(workflowDescriptor.workflowOptions.toGoogleAuthOptions), jesAttributes.endpointUrl)
  private def workflowPath: Path = jesCallPaths.workflowRootPath
  private def callRootPath: Path = jesCallPaths.callRootPath
  private def returnCodeFilename = jesCallPaths.returnCodeFilename
  private def returnCodeGcsPath = jesCallPaths.returnCodePath
  private def jesStdoutFile = jesCallPaths.stdoutPath
  private def jesStderrFile = jesCallPaths.stderrPath
  private def jesLogFilename = jesCallPaths.jesLogFilename
  private lazy val call = jobDescriptor.key.call
  private lazy val runtimeAttributes = {
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues evaluate
    // Fail the call if runtime attributes can't be evaluated
    val evaluatedAttributes = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    JesRuntimeAttributes(evaluatedAttributes)
  }
  override lazy val retryable = jobDescriptor.key.attempt <= runtimeAttributes.preemptible
  private lazy val workingDisk: JesAttachedDisk = runtimeAttributes.disks.find(_.name == JesWorkingDisk.Name).get
  private lazy val gcsExecPath: Path = callRootPath.resolve(JesExecScript)
  private lazy val cmdInput = JesFileInput(ExecParamName, gcsExecPath.toString, Paths.get(JesExecScript), workingDisk)
  private lazy val jesCommandLine = s"/bin/bash ${cmdInput.containerPath.toAbsolutePath.toString}"
  private lazy val defaultMonitoringOutputPath = callRootPath.resolve(JesMonitoringLogFile)
  private lazy val rcJesOutput = JesFileOutput(returnCodeFilename, returnCodeGcsPath.toString, Paths.get(returnCodeFilename), workingDisk)

  private lazy val standardParameters = Seq(rcJesOutput)
  // PBE this shouldn't need to be defined here, could be at the actor factory level.
  private lazy val jesAttributes = JesAttributes(googleConfig, configurationDescriptor)
  private lazy val returnCodeContents = Try(returnCodeGcsPath.toAbsolutePath.contentAsString)
  private lazy val dockerConfiguration: Option[JesDockerCredentials] = DockerConfiguration.build(backendConfig).dockerCredentials map JesDockerCredentials.apply
  private lazy val maxPreemption = runtimeAttributes.preemptible
  private lazy val preemptible: Boolean = jobDescriptor.key.attempt <= maxPreemption
  private lazy val tag = self.path.name

  private def globOutputPath(glob: String) = callRootPath.resolve(s"glob-${glob.md5Sum}/")

  private def gcsAuthFilePath: Path =  {
    // If we are going to upload an auth file we need a valid GCS path passed via workflow options.
    val bucket = workflowDescriptor.workflowOptions.get(WorkflowOptionKeys.AuthFilePath) getOrElse workflowPath.toString
    bucket.toPath(gcsFileSystem).resolve(s"${workflowDescriptor.id}_auth.json")
  }

  private def authGcsCredentialsPath(gcsPath: String): JesInput = JesLiteralInput(ExtraConfigParamName, gcsPath)

  private def gcsAuthParameter: Option[JesInput] = {
    if (jesAttributes.gcsFilesystemAuth.requiresAuthFile || dockerConfiguration.isDefined)
      Option(authGcsCredentialsPath(gcsAuthFilePath.toString))
    else None
  }

  private val callEngineFunctions = {
    val callContext = new CallContext(
      callRootPath,
      jesStdoutFile.toString,
      jesStderrFile.toString
    )

    new JesExpressionFunctions(List(gcsFileSystem), callContext)
  }

  private val lookup: ScopedLookupFunction = {
    val declarations = workflowDescriptor.workflowNamespace.workflow.declarations ++ call.task.declarations
    val unqualifiedWorkflowInputs = workflowDescriptor.inputs map {
      case (fqn, v) => splitFqn(fqn)._2 -> v
    }
    val knownInputs = unqualifiedWorkflowInputs ++ jobDescriptor.inputs
    WdlExpression.standardLookupFunction(knownInputs, declarations, callEngineFunctions)
  }

  private def evaluate(wdlExpression: WdlExpression) = wdlExpression.evaluate(lookup, callEngineFunctions)

  /**
    * Takes two arrays of remote and local WDL File paths and generates the necessary JesInputs.
    */
  private def jesInputsFromWdlFiles(jesNamePrefix: String,
                                    remotePathArray: Seq[WdlFile],
                                    localPathArray: Seq[WdlFile],
                                    jobDescriptor: BackendJobDescriptor): Iterable[JesInput] = {
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((remotePath, localPath), index) =>
        Seq(JesFileInput(s"$jesNamePrefix-$index", remotePath.valueString, Paths.get(localPath.valueString), workingDisk))
    }
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

  private def generateJesInputs(jobDescriptor: BackendJobDescriptor): Iterable[JesInput] = {
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
      val value = expression.evaluate(lookup, callEngineFunctions)
      expression.toWdlString.md5SumShort -> value
    } toMap

    val writeFunctionFiles = evaluatedExpressionMap collect { case (k, v: Success[_]) => k -> v.get } collect { case (k, v: WdlFile) => k -> Seq(v)}

    /** Collect all WdlFiles from inputs to the call */
    val callInputFiles: Map[FullyQualifiedName, Seq[WdlFile]] = jobDescriptor.inputs mapValues { _.collectAsSeq { case w: WdlFile => w } }

    (callInputFiles ++ writeFunctionFiles) flatMap {
      case (name, files) => jesInputsFromWdlFiles(name, files, files.map(relativeLocalizationPath), jobDescriptor)
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

  private def generateJesOutputs(jobDescriptor: BackendJobDescriptor): Seq[JesFileOutput] = {
    val wdlFileOutputs = jobDescriptor.key.scope.task.outputs flatMap { taskOutput =>
      taskOutput.requiredExpression.evaluateFiles(lookup, NoFunctions, taskOutput.wdlType) match {
        case Success(wdlFiles) => wdlFiles map relativeLocalizationPath
        case Failure(ex) =>
          log.warning(s"Could not evaluate $taskOutput: ${ex.getMessage}", ex)
          Seq.empty[WdlFile]
      }
    }

    // Create the mappings. GLOB mappings require special treatment (i.e. stick everything matching the glob in a folder)
    wdlFileOutputs.distinct map { wdlFile =>
      val destination = wdlFile match {
        case WdlSingleFile(filePath) => callRootPath.resolve(filePath).toString
        case WdlGlobFile(filePath) => globOutputPath(filePath).toString
      }
      val (relpath, disk) = relativePathAndAttachedDisk(wdlFile.value, runtimeAttributes.disks)
      JesFileOutput(makeSafeJesReferenceName(wdlFile.value), destination, relpath, disk)
    }
  }

  private def instantiateCommand: Try[String] = {
    val backendInputs = jobDescriptor.inputs mapValues gcsPathToLocal
    jobDescriptor.call.instantiateCommandLine(backendInputs, callEngineFunctions, gcsPathToLocal)
  }

  private def isFatalJesException(t: Throwable): Boolean = t match {
    case e: HttpResponseException if e.getStatusCode == 403 => true
    case _ => false
  }

  private def isTransientJesException(t: Throwable): Boolean = t match {
    // Quota exceeded
    case e: HttpResponseException if e.getStatusCode == 429 => true
    case _ => false
  }

  private def uploadCommandScript(command: String, withMonitoring: Boolean): Future[Unit] = {
    val monitoring = if (withMonitoring) {
      s"""|touch $JesMonitoringLogFile
          |chmod u+x $JesMonitoringScript
          |./$JesMonitoringScript > $JesMonitoringLogFile &""".stripMargin
    } else ""

    val tmpDir = Paths.get(JesWorkingDisk.MountPoint).resolve("tmp")
    val rcPath = Paths.get(JesWorkingDisk.MountPoint).resolve(returnCodeFilename)

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

    def writeScript(): Future[Unit] = Future(gcsExecPath.write(fileContent))

    implicit val system = context.system
    Retry.withRetry(
      writeScript,
      isTransient = isTransientJesException,
      isFatal = isFatalJesException
    )
  }

  private def googleProject(descriptor: BackendWorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleProject, jesAttributes.project)
  }

  private def createJesRun(jesParameters: Seq[JesParameter], runIdForResumption: Option[String] = None): Future[Run] = {

    def createRun() = Future(Run(
      runIdForResumption,
      jobDescriptor = jobDescriptor,
      runtimeAttributes = runtimeAttributes,
      callRootPath = callRootPath.toString,
      commandLine = jesCommandLine,
      logFileName = jesLogFilename,
      jesParameters,
      googleProject(jobDescriptor.descriptor),
      retryable,
      genomicsFactory
    ))

    implicit val system = context.system
    Retry.withRetry(
      createRun,
      isTransient = isTransientJesException,
      isFatal = isFatalJesException
    )
  }

  protected def runWithJes(command: String,
                           jesInputs: Seq[JesInput],
                           jesOutputs: Seq[JesFileOutput],
                           // runIdForResumption: Option[String],
                           withMonitoring: Boolean): Future[ExecutionHandle] = {

    val jesParameters = standardParameters ++ gcsAuthParameter ++ jesInputs ++ jesOutputs

    val jesJobSetup = for {
      _ <- uploadCommandScript(command, withMonitoring)
      run <- createJesRun(jesParameters)
    } yield run

    jesJobSetup map { run => JesPendingExecutionHandle(jobDescriptor, jesOutputs, run, previousStatus = None) }
  }

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    // PBE Currently assumes only execute
    val monitoringOutput = monitoringScript map { _ =>
      JesFileOutput(s"$MonitoringParamName-out", defaultMonitoringOutputPath.toString, Paths.get(JesMonitoringLogFile), workingDisk)
    }

    val jesInputs: Seq[JesInput] = generateJesInputs(jobDescriptor).toSeq ++ monitoringScript :+ cmdInput
    val jesOutputs: Seq[JesFileOutput] = generateJesOutputs(jobDescriptor) ++ monitoringOutput

    instantiateCommand match {
      case Success(command) =>
        runWithJes(command, jesInputs, jesOutputs, /* runIdForResumption, */ monitoringScript.isDefined)
      case Failure(ex: SocketTimeoutException) => Future.successful(FailedNonRetryableExecutionHandle(ex))
      case Failure(ex) => Future.successful(FailedNonRetryableExecutionHandle(ex))
    }
  }

  /**
    * Update the ExecutionHandle
    */
  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    previous match {
      case handle: JesPendingExecutionHandle =>
        val wfId = handle.jobDescriptor.descriptor.id.shortString
        val tag = handle.jobDescriptor.key.tag
        val runId = handle.run.runId
        log.debug(s"[UUID($wfId)$tag] Polling JES Job $runId")
        val status = Try(handle.run.checkStatus(jobDescriptor, handle.previousStatus))
        status match {
          case Success(s: TerminalRunStatus) => executionResult(s, handle)
          case Success(s) => handle.copy(previousStatus = Option(s)).future // Copy the current handle with updated previous status.
          case Failure(e: GoogleJsonResponseException) if e.getStatusCode == 404 =>
            log.error(s"JES Job ID ${handle.run.runId} has not been found, failing call")
            FailedNonRetryableExecutionHandle(e).future
          case Failure(e: Exception) =>
            // Log exceptions and return the original handle to try again.
            log.warning("Caught exception, retrying: " + e.getMessage, e)
            handle.future
          case Failure(e: Error) => Future.failed(e) // JVM-ending calamity.
          case Failure(throwable) =>
            // Someone has subclassed Throwable directly?
            FailedNonRetryableExecutionHandle(throwable).future
        }
      case f: FailedNonRetryableExecutionHandle => f.future
      case s: SuccessfulExecutionHandle => s.future
      case badHandle => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $badHandle"))
    }
  } flatten

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
  private def gcsInputToGcsOutput(inputValue: WdlValue): WdlValue = {
    // Convert to the local path where the file is localized to in the VM:
    val vmLocalizationPath = gcsPathToLocal(inputValue)

    vmLocalizationPath match {
      // If it's a file, work out where the file would be delocalized to, otherwise no-op:
      case x : WdlFile =>
        val delocalizationPath = callRootPath.resolve(vmLocalizationPath.valueString).toString
        WdlFile(delocalizationPath)
      case a: WdlArray => WdlArray(a.wdlType, a.value map { f => gcsInputToGcsOutput(f) })
      case m: WdlMap => WdlMap(m.wdlType, m.value map { case (k, v) => gcsInputToGcsOutput(k) -> gcsInputToGcsOutput(v) })
      case other => other
    }
  }

  private def customLookupFunction(alreadyGeneratedOutputs: Map[String, WdlValue]): String => WdlValue = toBeLookedUp => {
    val originalLookup: ScopedLookupFunction = x => alreadyGeneratedOutputs.getOrElse(x, lookup(x))
    gcsInputToGcsOutput(originalLookup(toBeLookedUp))
  }

  private def wdlValueToGcsPath(jesOutputs: Seq[JesFileOutput])(value: WdlValue): WdlValue = {
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

  private def outputLookup(taskOutput: TaskOutput, currentList: Seq[AttemptedLookupResult]) = for {
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
    wdlValue <- taskOutput.requiredExpression.evaluate(customLookupFunction(currentList.toLookupMap), callEngineFunctions)
    coercedValue <- taskOutput.wdlType.coerceRawValue(wdlValue)
    value = wdlValueToGcsPath(generateJesOutputs(jobDescriptor))(coercedValue)
  } yield value


  private def outputFoldingFunction: (Seq[AttemptedLookupResult], TaskOutput) => Seq[AttemptedLookupResult] = {
    (currentList: Seq[AttemptedLookupResult], taskOutput: TaskOutput) => {
      currentList ++ Seq(AttemptedLookupResult(taskOutput.name, outputLookup(taskOutput, currentList)))
    }
  }

  private def postProcess: Try[CallOutputs] = {
    val outputs = call.task.outputs
    val outputMappings = outputs.foldLeft(Seq.empty[AttemptedLookupResult])(outputFoldingFunction).map(_.toPair).toMap
    TryUtil.sequenceMap(outputMappings) map { outputMap =>
      outputMap mapValues { v =>
        // PBE fix hashing
        CallOutput(v, hash = None)
      }
    }
  }

  private def handleSuccess(outputMappings: Try[CallOutputs],
                            returnCode: Int,
                            hash: ExecutionHash,
                            executionHandle: ExecutionHandle): ExecutionHandle = {
    outputMappings match {
      case Success(outputs) => SuccessfulExecutionHandle(outputs, returnCode, hash)
      case Failure(ex: CromwellAggregatedException) if ex.throwables collectFirst { case s: SocketTimeoutException => s } isDefined =>
        // Return the execution handle in this case to retry the operation
        executionHandle
      case Failure(ex) => FailedNonRetryableExecutionHandle(ex)
    }
  }

  private def extractErrorCodeFromErrorMessage(errorMessage: String): Int = {
    errorMessage.substring(0, errorMessage.indexOf(':')).toInt
  }

  private def preempted(errorCode: Int, errorMessage: Option[String]): Boolean = {
    def isPreemptionCode(code: Int) = code == 13 || code == 14

    try {
      errorCode == 10 && errorMessage.isDefined && isPreemptionCode(extractErrorCodeFromErrorMessage(errorMessage.get)) && preemptible
    } catch {
      case _: NumberFormatException | _: StringIndexOutOfBoundsException =>
        log.warning(s"Unable to parse JES error code from error message: ${errorMessage.get}, assuming this was not a preempted VM.")
        false
    }
  }

  private def handleFailure(errorCode: Int, errorMessage: Option[String]) = {
    import lenthall.numeric.IntegerUtil._

    val taskName = s"${workflowDescriptor.id}:${call.unqualifiedName}"
    val attempt = jobDescriptor.key.attempt

    if (errorMessage.exists(_.contains("Operation canceled at")))  {
      AbortedExecutionHandle.future
    } else if (preempted(errorCode, errorMessage)) {
      val preemptedMsg = s"Task $taskName was preempted for the ${attempt.toOrdinal} time."

      if (attempt < maxPreemption) {
        val e = new PreemptedException(
          s"""$preemptedMsg The call will be restarted with another preemptible VM (max preemptible attempts number is $maxPreemption).
             |Error code $errorCode. Message: $errorMessage""".stripMargin
        )
        FailedRetryableExecutionHandle(e, None).future
      } else {
        val e = new PreemptedException(
          s"""$preemptedMsg The maximum number of preemptible attempts ($maxPreemption) has been reached. The call will be restarted with a non-preemptible VM.
             |Error code $errorCode. Message: $errorMessage)""".stripMargin)
        FailedRetryableExecutionHandle(e, None).future
      }
    } else {
      val e = new Throwable(s"Task ${workflowDescriptor.id}:${jobDescriptor.call.unqualifiedName} failed: error code $errorCode. Message: ${errorMessage.getOrElse("null")}")
      FailedNonRetryableExecutionHandle(e, None).future
    }
  }

  // PBE ideally hashes should be determinitic
  private def completelyRandomExecutionHash: Future[ExecutionHash] = Future.successful(ExecutionHash(UUID.randomUUID().toString, dockerHash = None))

  private def executionResult(status: RunStatus, handle: JesPendingExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    try {
      lazy val stderrLength: Long = jesStderrFile.size
      lazy val returnCode = returnCodeContents map { _.trim.toInt }
      lazy val continueOnReturnCode = runtimeAttributes.continueOnReturnCode

      status match {
        case RunStatus.Success if runtimeAttributes.failOnStderr && stderrLength.intValue > 0 =>
          // returnCode will be None if it couldn't be downloaded/parsed, which will yield a null in the DB
          FailedNonRetryableExecutionHandle(new Throwable(s"execution failed: stderr has length $stderrLength"), returnCode.toOption).future
        case RunStatus.Success if returnCodeContents.isFailure =>
          val exception = returnCode.failed.get
          log.warning(s"$tag could not download return code file, retrying: " + exception.getMessage, exception)
          // Return handle to try again.
          handle.future
        case RunStatus.Success if returnCode.isFailure =>
          FailedNonRetryableExecutionHandle(new Throwable(s"execution failed: could not parse return code as integer: " + returnCodeContents.get)).future
        case RunStatus.Success if !continueOnReturnCode.continueFor(returnCode.get) =>
          FailedNonRetryableExecutionHandle(new Throwable(s"execution failed: disallowed command return code: " + returnCode.get), returnCode.toOption).future
        case RunStatus.Success =>
          completelyRandomExecutionHash map { h => handleSuccess(postProcess, returnCode.get, h, handle) }
        case RunStatus.Failed(errorCode, errorMessage) => handleFailure(errorCode, errorMessage)
      }
    } catch {
      case e: Exception =>
        log.warning("Caught exception trying to download result, retrying: " + e.getMessage, e)
        // Return the original handle to try again.
        handle.future
    }
  } flatten

  protected implicit def ec: ExecutionContext = context.dispatcher
}
