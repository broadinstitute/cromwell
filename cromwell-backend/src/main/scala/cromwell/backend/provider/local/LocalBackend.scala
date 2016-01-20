package cromwell.backend.provider.local

import java.nio.file.{Files, Path, Paths}

import akka.actor.{Props, ActorRef}
import akka.util.Timeout
import better.files._
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.BackendActor
import cromwell.backend.model._
import cromwell.backend.provider.local.FileExtensions._
import wdl4s.types.WdlFileType
import wdl4s.values.WdlValue

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

  def props(task: TaskDescriptor): Props = Props(new LocalBackend(task))
}

/**
  * Executes a task in local computer through command line. It can be also executed in a Docker container.
  * @param task Task descriptor.
  */
class LocalBackend(task: TaskDescriptor) extends BackendActor with StrictLogging {

  import LocalBackend._

  implicit val timeout = Timeout(10 seconds)
  val subscriptions = ArrayBuffer[Subscription[ActorRef]]()

  val workingDir = task.workingDir
  val taskWorkingDir = task.name
  val executionDir = Paths.get(CromwellExecutionDir, workingDir, taskWorkingDir)
  val stdout = Paths.get(executionDir.toString, StdoutFile)
  val stderr = Paths.get(executionDir.toString, StderrFile)
  val script = Paths.get(executionDir.toString, ScriptFile)
  val returnCode = Paths.get(executionDir.toString, ReturnCodeFile)
  val stdoutWriter = stdout.untailed
  val stderrTailed = stderr.tailed(100)
  val argv = Seq("/bin/bash", script.toString)
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
    subscriptions.filter(subs => subs.eventType.isInstanceOf[ExecutionEvent]).foreach(
      subs => subs.subscriber ! new TaskStatus(Status.Created, None))
  }

  /**
    * Stops a task execution.
    */
  override def stop(): Unit = {
    processAbortFunc.get.apply()
    subscriptions.filter(subs => subs.eventType.isInstanceOf[ExecutionEvent]).foreach(
      subs => subs.subscriber ! new TaskStatus(Status.Canceled, None))
  }

  /**
    * Executes task in given context.
    */
  override def execute(): Unit = {
    subscriptions.filter(subs => subs.eventType.isInstanceOf[ExecutionEvent]).foreach(
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
    val sub = subscription.asInstanceOf[Subscription[ActorRef]]
    subscriptions += sub
    sub.subscriber ! Subscribed
  }

  /**
    * Unsubscribe to events on backend.
    */
  override def unsubscribeToEvent[T](subscription: Subscription[T]): Unit = {
    val sub = subscription.asInstanceOf[Subscription[ActorRef]]
    subscriptions -= sub
    sub.subscriber ! Unsubscribed
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
    val inputFolder = task.inputs.filter(p => p._2.wdlType.equals(WdlFileType)).flatMap(p => extractFolder(p._2.valueString)).toList.distinct
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
  private def executeTask(): TaskStatus = {
    def getCmdToExecute(): Seq[String] = dockerImage match {
      case Some(image) => buildDockerRunCommand(image).split(" ").toSeq
      case None => argv
    }

    val process = getCmdToExecute.run(ProcessLogger(stdoutWriter writeWithNewline, stderrTailed writeWithNewline))
    processAbortFunc = Some(() => process.destroy())
    subscriptions.filter(subs => subs.eventType.isInstanceOf[ExecutionEvent]).foreach(
      subs => subs.subscriber ! new TaskStatus(Status.Running, None))
    val backendCommandString = argv.map(s => "\"" + s + "\"").mkString(" ")
    logger.debug(s"command: $backendCommandString")
    val processReturnCode = process.exitValue() // blocks until process finishes

    Vector(stdoutWriter.writer, stderrTailed.writer) foreach {
      _.flushAndClose()
    }

    val stderrFileLength = Try(Files.size(stderr)).getOrElse(0L)

    if (stderrFileLength > 0) {
      // TODO: verify generated files exist
      TaskStatus(Status.Failed, Some(FailureResult(new IllegalStateException("StdErr file is not empty."), Some(processReturnCode), stderr.toString)))
    } else if (stderrFileLength <= 0 && processReturnCode != 0) {
      // TODO: verify generated files exist
      TaskStatus(Status.Failed, Some(FailureResult(new IllegalStateException("RC code is not equals to zero."), Some(processReturnCode), stderr.toString)))
    } else {
      val lookupFunction: String => WdlValue = inputName => task.inputs.get(inputName).get
      val outputsExpressions = task.outputs.map(output => output.name -> output.expression.evaluate(lookupFunction, new WorkflowEngineFunctions(executionDir)))
      if(outputsExpressions.filter(_._2.isFailure).size > 0)
        throw new IllegalStateException("Failed to evaluate output expressions!", outputsExpressions.filter(_._2.isFailure).head._2.failed.get)
      // TODO: verify generated files exist
      TaskStatus(Status.Succeeded, Some(SuccesfulResult(outputsExpressions.map(output => output._1 -> output._2.get).toMap)))
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