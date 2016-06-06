package cromwell.engine.backend.local

import java.io.Writer
import java.nio.file.{Files, Path, Paths}

import akka.actor.ActorSystem
import better.files._
import com.google.api.client.util.ExponentialBackOff.Builder
import cromwell.backend.JobKey
import cromwell.backend.wdl.OldCallEngineFunctions
import cromwell.core.PathFactory._
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.backend.local.OldStyleLocalBackend.InfoKeys
import org.slf4j.LoggerFactory
import wdl4s.values.WdlValue

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleLocalBackend {

  val ContainerRoot = "/root"
  val CallPrefix = "call"
  val ShardPrefix = "shard"

  object InfoKeys {
    val Pid = "PID"
  }

  /**
    * Simple utility implicit class adding a method that writes a line and appends a newline in one shot.
    */
  implicit class WriteWithNewline(val writer: Writer) extends AnyVal {
    def writeWithNewline(string: String): Unit = {
      writer.write(string)
      writer.write("\n")
    }
  }

  implicit class LocalEnhancedJobDescriptor(val jobDescriptor: OldStyleBackendCallJobDescriptor) extends AnyVal {
    def dockerContainerExecutionDir = jobDescriptor.workflowDescriptor.workflowRootPathWithBaseRoot(OldStyleLocalBackend.ContainerRoot)
    def containerCallRoot = jobDescriptor.callRuntimeAttributes.docker match {
      case Some(docker) => jobDescriptor.callRootPathWithBaseRoot(OldStyleLocalBackend.ContainerRoot)
      case None => jobDescriptor.callRootPath
    }
    def stdout = jobDescriptor.callRootPath.resolve("stdout")
    def stderr = jobDescriptor.callRootPath.resolve("stderr")
    def script = jobDescriptor.callRootPath.resolve("script")
    def returnCode = jobDescriptor.callRootPath.resolve("rc")
  }
}


/**
  * Handles both local Docker runs as well as local direct command line executions.
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class OldStyleLocalBackend(backendConfigEntry: BackendConfigurationEntry, actorSystem: ActorSystem) extends OldStyleBackend with SharedFileSystemBackend {

  private lazy val logger = LoggerFactory.getLogger("cromwell")

  override def adjustInputPaths(jobDescriptor: OldStyleBackendCallJobDescriptor) = adjustSharedInputPaths(jobDescriptor)

  /**
    * Exponential Backoff Builder to be used when polling for call status.
    */
  final private lazy val pollBackoffBuilder = new Builder()
    .setInitialIntervalMillis(10.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Int.MaxValue)
    .setMaxIntervalMillis(10.minutes.toMillis.toInt)
    .setMultiplier(1.1D)

  /**
    * If using Mac OS X with Docker Machine, by default Cromwell can only use the -v flag to 'docker run'
    * if running from somewhere within the user's home directory.  If not, then -v will mount to the
    * virtual machine that Docker Machine creates.  This ends up bubbling up to Cromwell when it can't find
    * an RC file in the call's output directory (because it's actually on the virtual machine).
    *
    * This function will at least detect this case and log some WARN messages.
    */
  def detectDockerMachinePossibleMisusage() = {
    val cromwellExecutionRoot = Paths.get(backendConfig.getString("root")).toAbsolutePath.toString
    val os = Option(System.getProperty("os.name"))
    val homeDir = Option(System.getProperty("user.home")).map(Paths.get(_).toAbsolutePath.toString)
    val dockerMachineName = sys.env.get("DOCKER_MACHINE_NAME")
    if (dockerMachineName.nonEmpty && os.contains("Mac OS X") && homeDir.isDefined && !cromwellExecutionRoot.startsWith(homeDir.get)) {
      val message = s"""You are running Docker using Docker Machine and your working directory is not a subdirectory of ${homeDir.get}
                        |By default, Docker Machine using Virtual Box only allows -v to mount to your Mac's filesystem for paths under your home directory.
                        |Running Cromwell outside of your home directory could lead to unexpected task failures.
                        |
                        |See https://docs.docker.com/engine/userguide/dockervolumes/ for more information on volume mounting in Docker."""
      message.stripMargin.split("\n").foreach(logger.warn)
    }
  }

  detectDockerMachinePossibleMisusage()

  override def pollBackoff = pollBackoffBuilder.build()

  def returnCode(jobDescriptor: OldStyleBackendCallJobDescriptor) = jobDescriptor.returnCode

  def stdoutStderr(jobDescriptor: OldStyleBackendCallJobDescriptor): CallLogs = sharedFileSystemStdoutStderr(jobDescriptor)

  def execute(jobDescriptor: OldStyleBackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = Future({
    val logger = jobLogger(jobDescriptor)
    jobDescriptor.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"`$instantiatedCommand`")
        writeScript(jobDescriptor, instantiatedCommand, jobDescriptor.containerCallRoot)
        runSubprocess(jobDescriptor)
      case Failure(ex) => OldStyleNonRetryableFailedExecution(ex).future
    }
  }).flatten map CompletedExecutionHandle

  override def isRestartable(key: JobKey, executionInfos: Map[String, Option[String]]): Boolean = true
  override def isResumable(key: JobKey, executionInfos: Map[String, Option[String]]): Boolean = false

  override def backendType = BackendType.LOCAL

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeScript(jobDescriptor: OldStyleBackendCallJobDescriptor, instantiatedCommand: String, containerRoot: Path) = {
    jobDescriptor.script.write(
      s"""#!/bin/sh
          |cd $containerRoot
          |$instantiatedCommand
          |echo $$? > rc
          |""".stripMargin)
  }

  /**
    * --rm automatically deletes the container upon exit
    * -v maps the host workflow executions directory to /root/<workflow id> on the container.
    * -i makes the run interactive, required for the cat and <&0 shenanigans that follow.
    */
  private def buildDockerRunCommand(jobDescriptor: OldStyleBackendCallJobDescriptor, image: String): String = {
    val callPath = jobDescriptor.containerCallRoot
    s"docker run --rm -v ${jobDescriptor.callRootPath.toAbsolutePath}:$callPath -i $image"
  }

  private def runSubprocess(jobDescriptor: OldStyleBackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[OldStyleExecutionResult] = {
    val logger = jobLogger(jobDescriptor)
    val stdoutWriter = jobDescriptor.stdout.untailed
    val stderrTailed = jobDescriptor.stderr.tailed(100)
    val dockerRun = jobDescriptor.callRuntimeAttributes.docker.map(d => buildDockerRunCommand(jobDescriptor, d)).getOrElse("")
    val argv = Seq("/bin/bash", "-c", s"cat ${jobDescriptor.script} | $dockerRun /bin/bash <&0")
    val process = argv.run(ProcessLogger(stdoutWriter writeWithNewline, stderrTailed writeWithNewline))

    // TODO: As currently implemented, this process.destroy() will kill the bash process but *not* its descendants. See ticket DSDEEPB-848.
    jobDescriptor.abortRegistrationFunction.foreach(_.register(AbortFunction(() => process.destroy())))
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"command: $backendCommandString")
    val processReturnCode = process.exitValue() // blocks until process finishes
    Vector(stdoutWriter.writer, stderrTailed.writer) foreach { _.flushAndClose() }

    val stderrFileLength = Try(Files.size(jobDescriptor.stderr)).getOrElse(0L)
    val returnCode = Try(
      if (processReturnCode == 0 || jobDescriptor.callRuntimeAttributes.docker.isEmpty) {
        val rc = jobDescriptor.returnCode.contentAsString.stripLineEnd.toInt
        logger.info(s"Return code: $rc")
        rc
      } else {
        logger.error(s"Non-zero return code: $processReturnCode")
        logger.error(s"Standard error was:\n\n${stderrTailed.tailString}\n")
        throw new Exception(s"Unexpected process exit code: $processReturnCode")
      }
    )

    if (jobDescriptor.callRuntimeAttributes.failOnStderr && stderrFileLength > 0) {
      OldStyleNonRetryableFailedExecution(new Throwable(s"Call ${jobDescriptor.call.fullyQualifiedName}, " +
        s"Workflow ${jobDescriptor.workflowDescriptor.id}: stderr has length $stderrFileLength")).future
    } else {

      def processSuccess(rc: Int) = {
        postProcess(jobDescriptor) match {
          case Success(outputs) => jobDescriptor.hash map { h => OldStyleSuccessfulBackendCallExecution(outputs, rc, h) }
          case Failure(e) =>
            val message = Option(e.getMessage) map { ": " + _ } getOrElse ""
            OldStyleNonRetryableFailedExecution(new Throwable("Failed post processing of outputs" + message, e)).future
        }
      }

      val badReturnCodeMessage =
        s"""Call ${jobDescriptor.call.fullyQualifiedName}, Workflow ${jobDescriptor.workflowDescriptor.id}: return code was ${returnCode.getOrElse("(none)")}
            |
            |Full command was: $backendCommandString
            |
            |${stderrTailed.tailString}
          """.stripMargin

      val continueOnReturnCode = jobDescriptor.callRuntimeAttributes.continueOnReturnCode
      returnCode match {
        case Success(143) => OldStyleAbortedExecution.future // Special case to check for SIGTERM exit code - implying abort
        case Success(otherReturnCode) if continueOnReturnCode.continueFor(otherReturnCode) => processSuccess(otherReturnCode)
        case Success(badReturnCode) => OldStyleNonRetryableFailedExecution(new Exception(badReturnCodeMessage), Option(badReturnCode)).future
        case Failure(e) => OldStyleNonRetryableFailedExecution(new Throwable(badReturnCodeMessage, e)).future
      }
    }
  }

  override def callRootPathWithBaseRoot(descriptor: OldStyleBackendCallJobDescriptor, baseRoot: String): Path = {
    val path = super.callRootPathWithBaseRoot(descriptor, baseRoot)
    if (!path.toFile.exists()) path.toFile.mkdirs()
    path
  }

  override def executionInfoKeys: List[String] = List(InfoKeys.Pid)

  override def callEngineFunctions(descriptor: OldStyleBackendCallJobDescriptor): OldCallEngineFunctions = {
    new OldStyleLocalCallEngineFunctions(descriptor.workflowDescriptor.fileSystems, buildCallContext(descriptor))
  }

  def instantiateCommand(jobDescriptor: OldStyleBackendCallJobDescriptor): Try[String] = {
    val backendInputs = adjustInputPaths(jobDescriptor)
    val pathTransformFunction: WdlValue => WdlValue = jobDescriptor.callRuntimeAttributes.docker match {
      case Some(_) => toDockerPath
      case None => (v: WdlValue) => v
    }
    jobDescriptor.call.instantiateCommandLine(backendInputs, jobDescriptor.callEngineFunctions, pathTransformFunction)
  }

  override def poll(jobDescriptor: OldStyleBackendCallJobDescriptor, previous: OldStyleExecutionHandle)(implicit ec: ExecutionContext) = Future.successful(previous)

  override def resume(descriptor: OldStyleBackendCallJobDescriptor, executionInfos: Map[String, Option[String]])(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = {
    Future.failed(new Throwable("resume invoked on non-resumable Local backend"))
  }
}
