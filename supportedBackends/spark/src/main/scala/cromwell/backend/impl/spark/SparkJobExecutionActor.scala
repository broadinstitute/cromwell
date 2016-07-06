package cromwell.backend.impl.spark

import java.nio.file.{FileSystems, Path}
import java.nio.file.attribute.PosixFilePermission

import akka.actor.Props
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, FailedNonRetryableResponse, SucceededResponse}
import cromwell.backend.io.{JobPaths, SharedFileSystem, SharedFsExpressionFunctions}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import wdl4s._
import wdl4s.types.WdlFileType
import wdl4s.util.TryUtil

import scala.concurrent.{Future, Promise}
import scala.sys.process.ProcessLogger
import scala.util.{Failure, Success, Try}
import scala.language.postfixOps

object SparkJobExecutionActor {
  val fileSystems = List(FileSystems.getDefault)

  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new SparkJobExecutionActor(jobDescriptor, configurationDescriptor))
}

class SparkJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                             override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor with SharedFileSystem {
  import SparkJobExecutionActor._
  import better.files._
  import cromwell.core.PathFactory._

  private val tag = s"CondorJobExecutionActor-${jobDescriptor.call.fullyQualifiedName}:"

  lazy val cmds = new SparkCommands
  lazy val extProcess = new SparkProcess
  private val fileSystemsConfig = configurationDescriptor.backendConfig.getConfig("filesystems")
  override val sharedFsConfig = fileSystemsConfig.getConfig("local")
  private val workflowDescriptor = jobDescriptor.descriptor
  private val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)

  // Files
  private val executionDir = jobPaths.callRoot
  private val scriptPath = jobPaths.script

  private lazy val stdoutWriter = extProcess.untailedWriter(jobPaths.stdout)
  private lazy val stderrWriter = extProcess.tailedWriter(100, jobPaths.stderr)

  private val call = jobDescriptor.key.call
  private val callEngineFunction = SharedFsExpressionFunctions(jobPaths, fileSystems)

  private val lookup = jobDescriptor.inputs.apply _

  private val executionResponse = Promise[BackendJobExecutionResponse]()

  private val runtimeAttributes = {
    val evaluateAttrs = call.task.runtimeAttributes.attrs mapValues (_.evaluate(lookup, callEngineFunction))
    // Fail the call if runtime attributes can't be evaluated
    val runtimeMap = TryUtil.sequenceMap(evaluateAttrs, "Runtime attributes evaluation").get
    SparkRuntimeAttributes(runtimeMap, jobDescriptor.descriptor.workflowOptions)
  }

  /**
    * Restart or resume a previously-started job.
    */
  override def recover: Future[BackendJobExecutionResponse] = {
    log.warning("{} HtCondor backend currently doesn't support recovering jobs. Starting {} again.", tag, jobDescriptor.key.call.fullyQualifiedName)
    Future(executeTask())
  }

  /**
    * Execute a new job.
    */
  override def execute: Future[BackendJobExecutionResponse] = {
    prepareAndExecute
    executionResponse.future
  }

  private def executeTask(): BackendJobExecutionResponse = {
    val argv = Seq(scriptPath.toString)
    val process = extProcess.externalProcess(argv, ProcessLogger(stdoutWriter.writeWithNewline, stderrWriter.writeWithNewline))
    val jobReturnCode = Try(process.exitValue()) // blocks until process (i.e. spark submission) finishes
    log.debug("{} Return code of spark submit command: {}", tag, jobReturnCode)
    List(stdoutWriter.writer, stderrWriter.writer).foreach(_.flushAndClose())
    jobReturnCode match {
      case Success(rc) if rc == 0 | runtimeAttributes.continueOnReturnCode.continueFor(rc) => processSuccess(rc)
      case Success(rc) => FailedNonRetryableResponse(jobDescriptor.key,
        new IllegalStateException(s"Execution process failed. Spark returned non zero status code: $jobReturnCode"), Option(rc))
      case Failure(error) => FailedNonRetryableResponse(jobDescriptor.key, error, None)
    }
  }

  private def processSuccess(rc: Int) = {
    evaluateOutputs(callEngineFunction, outputMapper(jobPaths)) match {
      case Success(outputs) => SucceededResponse(jobDescriptor.key, Some(rc), outputs)
      case Failure(e) =>
        val message = Option(e.getMessage) map {
          ": " + _
        } getOrElse ""
        FailedNonRetryableResponse(jobDescriptor.key, new Throwable("Failed post processing of outputs" + message, e), Option(rc))
    }
  }

  /**
    * Abort a running job.
    */
  override def abort(): Unit = Future.failed(new UnsupportedOperationException("SparkBackend currently doesn't support aborting jobs."))


  private def createExecutionFolderAndScript(): Unit = {
    try {
      log.debug("{} Creating execution folder: {}", tag, executionDir)
      executionDir.toString.toFile.createIfNotExists(true)

      log.debug("{} Resolving job command", tag)
      val command = localizeInputs(jobPaths.callRoot, docker = false, fileSystems, jobDescriptor.inputs) flatMap {
        localizedInputs => resolveJobCommand(localizedInputs)
      }

      log.debug("{} Creating bash script for executing command: {}", tag, command)
      cmds.writeScript(command.get, scriptPath, executionDir) // Writes the bash script for executing the command
      scriptPath.addPermission(PosixFilePermission.OWNER_EXECUTE) // Add executable permissions to the script.
    } catch {
      case ex: Exception =>
        log.error(ex, "Failed to prepare task: " + ex.getMessage)
        throw ex
    }
  }

  private def resolveJobCommand(localizedInputs: CallInputs): Try[String] = {
    if (runtimeAttributes.dockerImage.isDefined)
      modifyCommandForDocker(call.task.instantiateCommand(localizedInputs, callEngineFunction, identity), localizedInputs)
    else
      call.task.instantiateCommand(localizedInputs, callEngineFunction, identity)
  }

  private def modifyCommandForDocker(jobCmd: Try[String], localizedInputs: CallInputs): Try[String] = {
    Try {
      val inputFiles = localizedInputs.filter { case (k,v) => v.wdlType.equals(WdlFileType) }
      val dockerInputDataVol: Seq[String] = inputFiles.values.map { value =>
        val limit = value.valueString.lastIndexOf("/")
        value.valueString.substring(0, limit)
      } toSeq
      val dockerCmd = "docker run -w %s %s %s --rm %s %s"
      val dockerVolume = "-v %s:%s"
      val dockerVolumeInputs = s"$dockerVolume:ro"
      // `v.get` is safe below since we filtered the list earlier with only defined elements
      val inputVolumes = dockerInputDataVol.distinct.map(v => dockerVolumeInputs.format(v, v)).mkString(" ")
      val outputVolume = dockerVolume.format(executionDir.toAbsolutePath.toString, runtimeAttributes.dockerOutputDir.getOrElse(executionDir.toAbsolutePath.toString))
      val cmd = dockerCmd.format(runtimeAttributes.dockerWorkingDir.getOrElse(executionDir.toAbsolutePath.toString), inputVolumes, outputVolume, runtimeAttributes.dockerImage.get, jobCmd.get)
      log.debug(s"Docker command line to be used for task execution: $cmd.")
      cmd
    }
  }

  private def prepareAndExecute: Unit = {
    createExecutionFolderAndScript()
    executionResponse success executeTask()
  }


}
