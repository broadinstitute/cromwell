package cromwell.backend.impl.spark

import java.nio.file.{Path, FileSystems}
import java.nio.file.attribute.PosixFilePermission

import akka.actor.Props
import cromwell.backend.BackendJobExecutionActor.{SucceededResponse, FailedNonRetryableResponse, BackendJobExecutionResponse}
import cromwell.backend.io.{SharedFsExpressionFunctions, JobPaths, SharedFileSystem}
import cromwell.backend.{BackendJobExecutionActor, BackendConfigurationDescriptor, BackendJobDescriptor}
import wdl4s.util.TryUtil

import scala.concurrent.Future
import scala.sys.process.ProcessLogger
import scala.util.{Try, Failure, Success}

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

  lazy val extProcess = new SparkProcess
  private val fileSystemsConfig = configurationDescriptor.backendConfig.getConfig("filesystems")
  override val sharedFsConfig = fileSystemsConfig.getConfig("local")
  private val workflowDescriptor = jobDescriptor.descriptor
  private val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobDescriptor.key)

  // Files
  private val executionDir = jobPaths.callRoot
  private val scriptPath = jobPaths.script

  private val processStderr = executionDir.resolve("process.stderr")
  private val processStdout = executionDir.resolve("process.stdout")
  private lazy val stdoutWriter = extProcess.untailedWriter(processStdout)
  private lazy val stderrWriter = extProcess.tailedWriter(100, processStderr)

  private val call = jobDescriptor.key.call
  private val callEngineFunction = SharedFsExpressionFunctions(jobPaths, fileSystems)

  private val lookup = jobDescriptor.inputs.apply _

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
  override def execute: Future[BackendJobExecutionResponse] = Future(executeTask())

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
  override def abort(): Unit = Future.failed(new UnsupportedOperationException("HtCondorBackend currently doesn't support aborting jobs."))

  override def preStart(): Unit = {
    log.debug("{} Creating execution folder: {}", tag, executionDir)
    executionDir.toString.toFile.createIfNotExists(true)
    try {
      val command = localizeInputs(jobPaths.callRoot, docker = false, fileSystems, jobDescriptor.inputs) flatMap { localizedInputs =>
        call.task.instantiateCommand(localizedInputs, callEngineFunction, identity)
      }
      log.debug("{} Creating bash script for executing command: {}", tag, command)
      writeScript(command.get, scriptPath, executionDir) // Writes the bash script for executing the command
      scriptPath.addPermission(PosixFilePermission.OWNER_EXECUTE) // Add executable permissions to the script.
    } catch {
      case ex: Exception =>
        log.error(ex, "Failed to prepare task: " + ex.getMessage)
        throw ex
    }
  }

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeScript(instantiatedCommand: String, filePath: Path, containerRoot: Path) = {
    filePath.write(
      s"""#!/bin/sh
          |cd $containerRoot
          |$instantiatedCommand
          |echo $$? > rc
          |""".stripMargin)
  }
}
