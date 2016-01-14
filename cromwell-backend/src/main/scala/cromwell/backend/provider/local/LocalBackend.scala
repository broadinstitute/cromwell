package cromwell.backend.provider.local

import java.nio.file.{Files, Path, Paths}

import akka.actor.{ActorRef, ActorSystem}
import akka.util.Timeout
import better.files._
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.BackendActor
import cromwell.backend.model.OutputType.OutputType
import cromwell.backend.model._
import cromwell.backend.provider.local.FileExtensions._

import scala.collection.mutable.ArrayBuffer
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.sys.process._
import scala.util.Try

object LocalBackend {
  // Folders
  val CromwellExecutionDir = "cromwell-executions"

  // Files
  val ReturnCodeFile = "rc"
  val StdoutFile = "stdout"
  val StderrFile = "stderr"
  val ScriptFile = "script"

  val DockerFlag = "Docker"

  // don't know yet
  val CallPrefix = "call"
  val ShardPrefix = "shard"
}

/**
  * Executes a task in local computer through command line. It can be also executed in a Docker container.
  * @param task Task descriptor.
  * @param actorSystem Cromwell actor system.
  */
class LocalBackend(task: TaskDescriptor)(implicit actorSystem: ActorSystem) extends BackendActor with StrictLogging {

  import LocalBackend._

  implicit val timeout = Timeout(10 seconds)
  val subscriptions = ArrayBuffer[Subscription[ActorRef]]()

  val workingDir = task.workingDir
  val taskWorkingDir = workingDir + task.name
  val executionDir = Paths.get(CromwellExecutionDir, workingDir, taskWorkingDir)
  val stdout = Paths.get(executionDir.toString, StdoutFile)
  val stderr = Paths.get(executionDir.toString, StderrFile)
  val script = Paths.get(executionDir.toString, ScriptFile)
  val returnCode = Paths.get(executionDir.toString, ReturnCodeFile)
  val stdoutWriter = stdout.untailed
  val stderrTailed = stderr.tailed(100)
  val argv = Seq("/bin/bash", "-c", s"cat ${script}")
  val dockerImage = getDockerImage(task.runtimeAttributes)
  var processAbortFunc: Option[() => Unit] = None

  /**
    * Prepare the task and context for execution.
    */
  override def prepare(): Unit = {
    logger.debug(s"Creating execution folder: $executionDir")
    executionDir.toString.toFile.createIfNotExists(true)
    logger.debug(s"Creating bash script for executing command: ${task.command}.")
    writeBashScript(task.command, executionDir)
  }

  /**
    * Stops a task execution.
    */
  override def stop(): Unit = {
    processAbortFunc.get.apply()
  }

  /**
    * Executes task in given context.
    */
  override def execute(): Unit = {
    subscriptions.filter(subs => subs.eventType.equals(ExecutionEvent)).foreach(
      subs => subs.subscriber ! executeTask)
  }

  /**
    * Performs a cleanUp after the task was executed.
    */
  override def cleanUp(): Unit = {
    // Do nothing.
  }

  /**
    * Subscribe to events on backend.
    */
  override def subscribeToEvent[T](subscription: Subscription[T]): Unit = {
    subscriptions += subscription.asInstanceOf[Subscription[ActorRef]]
  }

  /**
    * Unsubscribe to events on backend.
    */
  override def unsubscribeToEvent[T](subscription: Subscription[T]): Unit = {
    subscriptions -= subscription.asInstanceOf[Subscription[ActorRef]]
  }

  /**
    * Get file content hash.
    * @param files List of files to compute.
    * @return List of hashes related to files.
    */
  override def computeInputFileHash(files: List[Path]): Map[Path, String] = {
    files.map(path => path -> File(path.toString).md5).toMap[Path, String]
  }

  /**
    * Get container image hash.
    * @param imageName Image name.
    * @return A hash.
    */
  override def computeContainerImageHash(imageName: String): Map[String, String] = ???

  /**
    * Gather Docker image name from runtime attributes. It it's not present returns none.
    * @param runtimeAttributes Runtime requirements for the task.
    */
  private def getDockerImage(runtimeAttributes: Map[String, String]): Option[String] = {
    if (runtimeAttributes.contains(DockerFlag) && runtimeAttributes.filter(p => p._1.equals(DockerFlag)).size == 1 &&
      !runtimeAttributes.filter(p => p._1.equals(DockerFlag)).values.head.isEmpty) {
      Some(runtimeAttributes.filter(p => p._1.equals(DockerFlag)).values.head)
    } else {
      None
    }
  }

  /**
    * Creates docker command in order to execute the task into a container.
    * @param image Docker image name.
    * @return Command to execute.
    */
  private def buildDockerRunCommand(image: String): String = {
    def extractFolder(file: String): String = file.substring(0, file.lastIndexOf("/"))
    val dockerVolume = "-v %s:%s"
    val inputFolder = task.inputs.filter(p => p._1.equals(OutputType.File)).flatMap(p => extractFolder(p._2)).toList.distinct
    val inputVolumes = inputFolder.map(v => dockerVolume.format(v, v)).mkString(" ")
    val outputVolume = dockerVolume.format(executionDir, executionDir)
    log.debug(s"DockerInputVolume: $inputVolumes")
    log.debug(s"DockerOutputVolume: $outputVolume")
    val dockerCmd = s"docker run %s %s --rm %s %s" //TODO: make it configurable.
    dockerCmd.format(inputVolumes, outputVolume, image, argv.mkString(" "))
  }

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  private def writeBashScript(taskCommand: String, executionPath: Path): Unit = {
    script.write(
      s"""#!/bin/sh
          |cd $executionPath
          |$taskCommand
          |echo $$? > rc
          |""".stripMargin)
  }

  /**
    * Run a command using a bash script.
    * @return A TaskStatus with the final status of the task.
    */
  private def executeTask(): TaskStatus[Option[Map[OutputType, String]]] = {
    def getCmdToExecute(): Seq[String] = dockerImage match {
      case Some(image) => buildDockerRunCommand(image).split(" ").toSeq
      case None => argv
    }

    val process = getCmdToExecute.run(ProcessLogger(stdoutWriter writeWithNewline, stderrTailed writeWithNewline))
    processAbortFunc = Some(() => process.destroy())

    val backendCommandString = argv.map(s => "\"" + s + "\"").mkString(" ")
    logger.debug(s"command: $backendCommandString")
    val processReturnCode = process.exitValue() // blocks until process finishes

    Vector(stdoutWriter.writer, stderrTailed.writer) foreach {
      _.flushAndClose()
    }

    val stderrFileLength = Try(Files.size(stderr)).getOrElse(0L)

    if (stderrFileLength > 0) {
      // TODO: verify generated files exist
      TaskStatus(Status.Failed, Option(Map(OutputType.StdError -> stderr.toString)))
    } else if (stderrFileLength <= 0 && processReturnCode != 0) {
      // TODO: verify generated files exist
      TaskStatus(Status.Failed, Option(Map(OutputType.StdError -> returnCode.toString)))
    } else {
      // TODO: verify generated files exist
      TaskStatus(Status.Succeeded, Option(task.outputs))
    }
  }

  // TODO: check continue on error.
  /*
    val badReturnCodeMessage =
      s"""Call ${backendCall.call.fullyQualifiedName}, Workflow ${backendCall.workflowDescriptor.id}: return code was ${returnCode.getOrElse("(none)")}
          |
          |Full command was: $backendCommandString
          |
          |${stderrTailed.tailString}
        """.stripMargin

    val continueOnReturnCode = backendCall.runtimeAttributes.continueOnReturnCode
    returnCode match {
      case Success(143) => AbortedExecution.future // Special case to check for SIGTERM exit code - implying abort
      case Success(otherReturnCode) if continueOnReturnCode.continueFor(otherReturnCode) => processSuccess(otherReturnCode)
      case Success(badReturnCode) => FailedExecution(new Exception(badReturnCodeMessage), Option(badReturnCode)).future
      case Failure(e) => FailedExecution(new Throwable(badReturnCodeMessage, e)).future
    }
  }*/
}