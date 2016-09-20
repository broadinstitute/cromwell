package cromwell.backend.impl.jes

import java.net.SocketTimeoutException
import java.nio.file.{Path, Paths}
import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.event.LoggingReceive
import better.files._
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.{ExecutionMode, JobId}
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, FailedRetryableExecutionHandle, NonRetryableExecution, SuccessfulExecutionHandle}
import cromwell.backend.impl.jes.JesImplicits.PathString
import cromwell.backend.impl.jes.JesJobExecutionActor.JesOperationIdKey
import cromwell.backend.impl.jes.RunStatus.TerminalRunStatus
import cromwell.backend.impl.jes.io._
import cromwell.backend.{AttemptedLookupResult, BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor, PreemptedException}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core._
import cromwell.core.logging.JobLogging
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}
import cromwell.filesystems.gcs.NioGcsPath
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.services.metadata.CallMetadataKeys._
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata._
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

  def props(jobDescriptor: BackendJobDescriptor,
            completionPromise: Promise[BackendJobExecutionResponse],
            jesWorkflowInfo: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef): Props = {
    Props(new JesAsyncBackendJobExecutionActor(jobDescriptor,
      completionPromise,
      jesWorkflowInfo,
      initializationData,
      serviceRegistryActor)).withDispatcher(BackendDispatcher)
  }

  object WorkflowOptionKeys {
    val MonitoringScript = "monitoring_script"
    val GoogleProject = "google_project"
  }

  private val ExecParamName = "exec"
  private val MonitoringParamName = "monitoring"

  private val JesMonitoringScript = JesWorkingDisk.MountPoint.resolve("monitoring.sh")
  private val JesMonitoringLogFile = JesWorkingDisk.MountPoint.resolve("monitoring.log")
  private val JesExecScript = "exec.sh"
  private val ExtraConfigParamName = "__extra_config_gcs_path"

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

  case class JesJobId(operationId: String) extends JobId
}


class JesAsyncBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                       override val completionPromise: Promise[BackendJobExecutionResponse],
                                       jesConfiguration: JesConfiguration,
                                       initializationData: JesBackendInitializationData,
                                       serviceRegistryActor: ActorRef)
  extends Actor with ActorLogging with AsyncBackendJobExecutionActor with JobLogging {

  import JesAsyncBackendJobExecutionActor._

  override lazy val pollBackOff = SimpleExponentialBackoff(
    initialInterval = 30 seconds, maxInterval = 10 minutes, multiplier = 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(
    initialInterval = 3 seconds, maxInterval = 20 seconds, multiplier = 1.1)

  private lazy val workflowDescriptor = jobDescriptor.workflowDescriptor

  // For Logging
  val workflowId = workflowDescriptor.id
  val jobTag = jobDescriptor.key.tag

  private lazy val jesAttributes = jesConfiguration.jesAttributes
  private[jes] lazy val jesCallPaths = initializationData.workflowPaths.toJesCallPaths(jobDescriptor.key)
  private[jes] lazy val monitoringScript: Option[JesInput] = {
    jobDescriptor.workflowDescriptor.workflowOptions.get(WorkflowOptionKeys.MonitoringScript) map { path =>
      JesFileInput(s"$MonitoringParamName-in", getPath(path).toString, JesWorkingDisk.MountPoint.resolve(JesMonitoringScript), workingDisk)
    } toOption
  }

  private def callRootPath: Path = jesCallPaths.callRootPath
  private def returnCodeFilename = jesCallPaths.returnCodeFilename
  private def returnCodeGcsPath = jesCallPaths.returnCodePath
  private def jesStdoutFile = jesCallPaths.stdoutPath
  private def jesStderrFile = jesCallPaths.stderrPath
  private def jesLogFilename = jesCallPaths.jesLogFilename
  private lazy val call = jobDescriptor.key.call
  private lazy val runtimeAttributes = JesRuntimeAttributes(jobDescriptor.runtimeAttributes, jobLogger)

  override lazy val retryable = jobDescriptor.key.attempt <= runtimeAttributes.preemptible
  private lazy val workingDisk: JesAttachedDisk = runtimeAttributes.disks.find(_.name == JesWorkingDisk.Name).get
  private lazy val gcsExecPath: Path = callRootPath.resolve(JesExecScript)
  private lazy val cmdInput = JesFileInput(ExecParamName, gcsExecPath.toString, Paths.get(JesExecScript), workingDisk)
  private lazy val jesCommandLine = s"/bin/bash ${cmdInput.containerPath}"
  private lazy val defaultMonitoringOutputPath = callRootPath.resolve(JesMonitoringLogFile)
  private lazy val rcJesOutput = JesFileOutput(returnCodeFilename, returnCodeGcsPath.toString, Paths.get(returnCodeFilename), workingDisk)

  private lazy val standardParameters = Seq(rcJesOutput)
  private lazy val returnCodeContents = Try(File(returnCodeGcsPath).contentAsString)
  private lazy val dockerConfiguration = jesConfiguration.dockerCredentials
  private lazy val maxPreemption = runtimeAttributes.preemptible
  private[jes] lazy val preemptible: Boolean = jobDescriptor.key.attempt <= maxPreemption
  private lazy val tag = s"${this.getClass.getSimpleName} [UUID(${workflowId.shortString}):${jobDescriptor.key.tag}]"

  private var runId: Option[String] = None

  private lazy val metadataJobKey = {
    val jobDescriptorKey: BackendJobDescriptorKey = jobDescriptor.key
    MetadataJobKey(jobDescriptorKey.call.fullyQualifiedName, jobDescriptorKey.index, jobDescriptorKey.attempt)
  }

  def jesReceiveBehavior: Receive = LoggingReceive {
    case AbortJobCommand =>
      runId foreach { id =>
        Try(Run(id, initializationData.genomics).abort()) match {
          case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, id)
          case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, id, ex.getMessage)
        }
      }
      context.parent ! AbortedResponse(jobDescriptor.key)
      context.stop(self)
    case KvPutSuccess(_) => // expected after the KvPut for the operation ID
  }

  override def receive: Receive = jesReceiveBehavior orElse super.receive

  private def globOutputPath(glob: String) = callRootPath.resolve(s"glob-${glob.md5Sum}/")

  private def gcsAuthParameter: Option[JesInput] = {
    if (jesAttributes.gcsFilesystemAuth.requiresAuthFile || dockerConfiguration.isDefined)
      Option(JesLiteralInput(ExtraConfigParamName, jesCallPaths.gcsAuthFilePath.toString))
    else None
  }

  private lazy val callContext = CallContext(
    callRootPath,
    jesStdoutFile.toString,
    jesStderrFile.toString
  )

  private[jes] lazy val callEngineFunctions = new JesExpressionFunctions(List(jesCallPaths.gcsFileSystem), callContext)

  private val lookup: ScopedLookupFunction = {
    val declarations = workflowDescriptor.workflowNamespace.workflow.declarations ++ call.task.declarations
    WdlExpression.standardLookupFunction(jobDescriptor.inputs, declarations, callEngineFunctions)
  }

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
    Try(getPath(file.value)) match {
      case Success(gcsPath: NioGcsPath) => WdlFile(gcsPath.bucket + "/" + gcsPath.objectName, file.isGlob)
      case Success(gcsPath) => file
      case Failure(e) => file
    }
  }

  private[jes] def generateJesInputs(jobDescriptor: BackendJobDescriptor): Iterable[JesInput] = {
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
      case p if !p.isAbsolute => JesWorkingDisk.MountPoint.resolve(p)
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

  private[jes] def generateJesOutputs(jobDescriptor: BackendJobDescriptor): Seq[JesFileOutput] = {
    val wdlFileOutputs = jobDescriptor.key.scope.task.outputs flatMap { taskOutput =>
      taskOutput.requiredExpression.evaluateFiles(lookup, NoFunctions, taskOutput.wdlType) match {
        case Success(wdlFiles) => wdlFiles map relativeLocalizationPath
        case Failure(ex) =>
          jobLogger.warn(s"Could not evaluate $taskOutput: ${ex.getMessage}", ex)
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

  private def uploadCommandScript(command: String, withMonitoring: Boolean): Future[Unit] = {
    val monitoring = if (withMonitoring) {
      s"""|touch $JesMonitoringLogFile
          |chmod u+x $JesMonitoringScript
          |$JesMonitoringScript > $JesMonitoringLogFile &""".stripMargin
    } else ""

    val tmpDir = File(JesWorkingDisk.MountPoint)./("tmp").path
    val rcPath = File(JesWorkingDisk.MountPoint)./(returnCodeFilename).path

    val fileContent =
      s"""
         |#!/bin/bash
         |export _JAVA_OPTIONS=-Djava.io.tmpdir=$tmpDir
         |export TMPDIR=$tmpDir
         |$monitoring
         |(
         |cd ${JesWorkingDisk.MountPoint}
         |$command
         |)
         |echo $$? > $rcPath
       """.stripMargin.trim

    def writeScript(): Future[Unit] = Future(File(gcsExecPath).write(fileContent))

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
      googleProject(jobDescriptor.workflowDescriptor),
      retryable,
      initializationData.genomics
    ))

    implicit val system = context.system
    Retry.withRetry(
      createRun,
      isTransient = isTransientJesException,
      isFatal = isFatalJesException
    ) andThen {
      case Success(run) =>
        // If this execution represents a resumption don't publish the operation ID since clearly it is already persisted.
        runId = Option(run.runId)
        if (runIdForResumption.isEmpty) {
          serviceRegistryActor ! KvPut(KvPair(ScopedKey(jobDescriptor.workflowDescriptor.id,
            KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt),
            JesOperationIdKey), runId))
        }
    }
  }

  protected def runWithJes(command: String,
                           jesInputs: Seq[JesInput],
                           jesOutputs: Seq[JesFileOutput],
                           runIdForResumption: Option[String],
                           withMonitoring: Boolean): Future[ExecutionHandle] = {

    tellStartMetadata()

    val jesParameters = standardParameters ++ gcsAuthParameter ++ jesInputs ++ jesOutputs

    val jesJobSetup = for {
      _ <- uploadCommandScript(command, withMonitoring)
      run <- createJesRun(jesParameters, runIdForResumption)
      _ = tellRunMetadata(run)
    } yield run

    jesJobSetup map { run => JesPendingExecutionHandle(jobDescriptor, jesOutputs, run, previousStatus = None) }
  }

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    val monitoringOutput = monitoringScript map { _ => JesFileOutput(s"$MonitoringParamName-out",
      defaultMonitoringOutputPath.toString, File(JesMonitoringLogFile).path, workingDisk)
    }

    tellMetadata(CallMetadataKeys.CallRoot, callRootPath)
    tellMetadata(JesMetadataKeys.GoogleProject, jesAttributes.project)
    tellMetadata(JesMetadataKeys.ExecutionBucket, jesAttributes.executionBucket)
    tellMetadata(JesMetadataKeys.EndpointUrl, jesAttributes.endpointUrl)
    monitoringOutput foreach { o => tellMetadata(JesMetadataKeys.MonitoringLog, o.gcs) }

    // Force runtimeAttributes to evaluate so we can fail quickly now if we need to:
    Try(runtimeAttributes) match {
      case Success(_) => startExecuting(monitoringOutput, mode)
      case Failure(e) => Future.successful(FailedNonRetryableExecutionHandle(e, None))
    }
  }

  private def startExecuting(monitoringOutput: Option[JesFileOutput], mode: ExecutionMode): Future[ExecutionHandle] = {
    val jesInputs: Seq[JesInput] = generateJesInputs(jobDescriptor).toSeq ++ monitoringScript :+ cmdInput
    val jesOutputs: Seq[JesFileOutput] = generateJesOutputs(jobDescriptor) ++ monitoringOutput

    instantiateCommand match {
      case Success(command) => runWithJes(command, jesInputs, jesOutputs, mode.jobId.collectFirst { case j: JesJobId => j.operationId }, monitoringScript.isDefined)
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
        val runId = handle.run.runId
        jobLogger.debug(s"$tag Polling JES Job $runId")
        val previousStatus = handle.previousStatus
        val status = Try(handle.run.status())
        status foreach { currentStatus =>
          if (!(handle.previousStatus contains currentStatus)) {
            // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise
            // just use the state names.
            val prevStateName = previousStatus map { _.toString } getOrElse "-"
            jobLogger.info(s"$tag Status change from $prevStateName to $currentStatus")
            tellMetadata("backendStatus", currentStatus.toString)
          }
        }
        status match {
          case Success(s: TerminalRunStatus) =>
            tellEventMetadata(s.eventList)
            tellMetadata(JesMetadataKeys.MachineType, s.machineType.getOrElse("unknown"))
            tellMetadata(JesMetadataKeys.InstanceName, s.instanceName.getOrElse("unknown"))
            tellMetadata(JesMetadataKeys.Zone, s.zone.getOrElse("unknown"))
            executionResult(s, handle)
          case Success(s) => handle.copy(previousStatus = Option(s)).future // Copy the current handle with updated previous status.
          case Failure(e: GoogleJsonResponseException) if e.getStatusCode == 404 =>
            jobLogger.error(s"$tag JES Job ID ${handle.run.runId} has not been found, failing call")
            FailedNonRetryableExecutionHandle(e).future
          case Failure(e: Exception) =>
            // Log exceptions and return the original handle to try again.
            jobLogger.warn(s"Caught exception, retrying", e)
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
    * Fire and forget start info to the metadata service
    */
  private def tellStartMetadata(): Unit = {
    val runtimeAttributesEvent = runtimeAttributes.asMap map {
      case (key, value) => MetadataEvent(metadataKey(s"runtimeAttributes:$key"), MetadataValue(value))
    }

    val events = runtimeAttributesEvent ++ List(
      metadataEvent("preemptible", preemptible),
      metadataEvent(Stdout, jesCallPaths.stdoutPath),
      metadataEvent(Stderr, jesCallPaths.stderrPath),
      metadataEvent(BackendLogsPrefix + ":log", jesCallPaths.jesLogPath),
      metadataEvent("cache:allowResultReuse", true)
    )

    serviceRegistryActor ! PutMetadataAction(events)
  }

  /**
    * Fire and forget events to the metadata service
    */
  private def tellEventMetadata(eventList: Seq[EventStartTime]): Unit = {
    eventList.headOption foreach { firstEvent =>
      // The final event is only used as the book-end for the final pairing so the name is never actually used...
      val offset = firstEvent.offsetDateTime.getOffset
      val now = OffsetDateTime.now.withOffsetSameInstant(offset)
      val lastEvent = EventStartTime("unused_name", now)
      val tailedEventList = eventList :+ lastEvent
      val events = tailedEventList.sliding(2).zipWithIndex flatMap {
        case (Seq(eventCurrent, eventNext), index) =>
          val eventKey = s"executionEvents[$index]"
          List(
            metadataEvent(s"$eventKey:description", eventCurrent.name),
            metadataEvent(s"$eventKey:startTime", eventCurrent.offsetDateTime),
            metadataEvent(s"$eventKey:endTime", eventNext.offsetDateTime)
          )
      }

      serviceRegistryActor ! PutMetadataAction(events.toIterable)
    }
  }

  /**
    * Fire and forget run info to the metadata service
    */
  private def tellRunMetadata(run: Run): Unit = {
    tellMetadata("jobId", run.runId)
  }

  /**
    * Fire and forget to the metadata service
    */
  private def tellMetadata(key: String, value: Any): Unit = {
    val event = metadataEvent(key, value)
    val putMetadataAction = PutMetadataAction(event)
    serviceRegistryActor ! putMetadataAction // and forget
  }

  private def metadataEvent(key: String, value: Any) = {
    val metadataValue = MetadataValue(value)
    MetadataEvent(metadataKey(key), metadataValue)
  }

  private def metadataKey(key: String) = MetadataKey(workflowId, Option(metadataJobKey), key)

  private def customLookupFunction(alreadyGeneratedOutputs: Map[String, WdlValue])(toBeLookedUp: String): WdlValue = alreadyGeneratedOutputs.getOrElse(toBeLookedUp, lookup(toBeLookedUp))

  private[jes] def wdlValueToGcsPath(jesOutputs: Seq[JesFileOutput])(value: WdlValue): WdlValue = {
    def toGcsPath(wdlFile: WdlFile) = jesOutputs collectFirst {
      case o if o.name == makeSafeJesReferenceName(wdlFile.valueString) => WdlFile(o.gcs)
    } getOrElse value
    value match {
      case wdlArray: WdlArray => wdlArray map wdlValueToGcsPath(jesOutputs)
      case wdlMap: WdlMap => wdlMap map {
        case (k, v) => wdlValueToGcsPath(jesOutputs)(k) -> wdlValueToGcsPath(jesOutputs)(v)
      }
      case file: WdlFile => if (file.value.isGcsUrl) file else toGcsPath(file)
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

  private def postProcess: Try[JobOutputs] = {
    val outputs = call.task.outputs
    val outputMappings = outputs.foldLeft(Seq.empty[AttemptedLookupResult])(outputFoldingFunction).map(_.toPair).toMap
    TryUtil.sequenceMap(outputMappings) map { outputMap =>
      outputMap mapValues { v => JobOutput(v) }
    }
  }

  private def gatherJobDetritusFiles: Map[String,String] = {
      Map("stdout" -> jesCallPaths.stdoutPath.toString, "stderr" -> jesCallPaths.stderrPath.toString,
        "returnCode" -> jesCallPaths.returnCodePath.toString, "jesLog" -> jesCallPaths.jesLogPath.toString,
        "gcsExec" -> gcsExecPath.toString, jesCallPaths.CallRootPathKey -> callRootPath.toString)
  }

  private def handleSuccess(outputMappings: Try[JobOutputs], returnCode: Int, jobDetritusFiles: Map[String, String], executionHandle: ExecutionHandle): ExecutionHandle = {
    outputMappings match {
      case Success(outputs) => SuccessfulExecutionHandle(outputs, returnCode, jobDetritusFiles)
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
        jobLogger.warn(s"Unable to parse JES error code from error message: {}, assuming this was not a preempted VM.", errorMessage.get)
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
        val e = PreemptedException(
          s"""$preemptedMsg The call will be restarted with another preemptible VM (max preemptible attempts number is $maxPreemption).
             |Error code $errorCode. Message: $errorMessage""".stripMargin
        )
        FailedRetryableExecutionHandle(e, None).future
      } else {
        val e = PreemptedException(
          s"""$preemptedMsg The maximum number of preemptible attempts ($maxPreemption) has been reached. The call will be restarted with a non-preemptible VM.
             |Error code $errorCode. Message: $errorMessage)""".stripMargin)
        FailedRetryableExecutionHandle(e, None).future
      }
    } else {
      val id = workflowDescriptor.id
      val name = jobDescriptor.call.unqualifiedName
      val message = errorMessage.getOrElse("null")
      val exception = new RuntimeException(s"Task $id:$name failed: error code $errorCode. Message: $message")
      FailedNonRetryableExecutionHandle(exception, None).future
    }
  }

  private[jes] def executionResult(status: RunStatus, handle: JesPendingExecutionHandle)
                                  (implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    try {
      lazy val stderrLength: Long = File(jesStderrFile).size
      lazy val returnCode = returnCodeContents map { _.trim.toInt }
      lazy val continueOnReturnCode = runtimeAttributes.continueOnReturnCode

      status match {
        case _: RunStatus.Success if runtimeAttributes.failOnStderr && stderrLength.intValue > 0 =>
          // returnCode will be None if it couldn't be downloaded/parsed, which will yield a null in the DB
          FailedNonRetryableExecutionHandle(new RuntimeException(
            s"execution failed: stderr has length $stderrLength"), returnCode.toOption).future
        case _: RunStatus.Success if returnCodeContents.isFailure =>
          val exception = returnCode.failed.get
          jobLogger.warn(s"could not download return code file, retrying", exception)
          // Return handle to try again.
          handle.future
        case _: RunStatus.Success if returnCode.isFailure =>
          FailedNonRetryableExecutionHandle(new RuntimeException(
            s"execution failed: could not parse return code as integer: ${returnCodeContents.get}")).future
        case _: RunStatus.Success if !continueOnReturnCode.continueFor(returnCode.get) =>
          val badReturnCodeMessage = s"Call ${call.fullyQualifiedName}: return code was ${returnCode.getOrElse("(none)")}"
          FailedNonRetryableExecutionHandle(new RuntimeException(badReturnCodeMessage), returnCode.toOption).future
        case _: RunStatus.Success => handleSuccess(postProcess, returnCode.get, gatherJobDetritusFiles, handle).future
        case RunStatus.Failed(errorCode, errorMessage, _, _, _, _) => handleFailure(errorCode, errorMessage)
      }
    } catch {
      case e: Exception =>
        jobLogger.warn("Caught exception trying to download result, retrying", e)
        // Return the original handle to try again.
        handle.future
    }
  } flatten

  /**
    * Takes a path in GCS and comes up with a local path which is unique for the given GCS path.
    *
    * Matches the path generated via relativeLocalizationPath and passed in as JesFileInput.local.
    *
    * @param mountPoint The mount point for inputs
    * @param gcsPath The input path
    * @return A path which is unique per input path
    */
  private def localFilePathFromCloudStoragePath(mountPoint: Path, gcsPath: NioGcsPath): Path = {
    mountPoint.resolve(gcsPath.bucket).resolve(gcsPath.objectName)
  }

  /**
    * Takes a single WdlValue and maps google cloud storage (GCS) paths into an appropriate local file path.
    * If the input is not a WdlFile, or the WdlFile is not a GCS path, the mapping is a noop.
    *
    * @param wdlValue the value of the input
    * @return a new FQN to WdlValue pair, with WdlFile paths modified if appropriate.
    */
  private[jes] def gcsPathToLocal(wdlValue: WdlValue): WdlValue = {
    wdlValue match {
      case wdlFile: WdlFile =>
        Try(getPath(wdlFile.valueString)) match {
          case Success(gcsPath: NioGcsPath) =>
            WdlFile(localFilePathFromCloudStoragePath(workingDisk.mountPoint, gcsPath).toString, wdlFile.isGlob)
          case Success(otherPath) => wdlValue
          case Failure(e) => wdlValue
        }
      case wdlArray: WdlArray => wdlArray map gcsPathToLocal
      case wdlMap: WdlMap => wdlMap map { case (k, v) => gcsPathToLocal(k) -> gcsPathToLocal(v) }
      case _ => wdlValue
    }
  }

  private def getPath(str: String) = jesCallPaths.gcsFileSystem.getPath(str)

  protected implicit def ec: ExecutionContext = context.dispatcher
}
