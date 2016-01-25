package cromwell.backend.provider.local

import java.nio.file.{Files, Path, Paths}
import java.util.regex.Pattern

import akka.actor.{ActorRef, Props}
import akka.util.Timeout
import better.files._
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.BackendActor
import cromwell.backend.model._
import cromwell.backend.provider.local.FileExtensions._
import wdl4s.WdlExpression
import wdl4s.types.{WdlFileType, WdlType}
import wdl4s.values.{WdlFile, WdlSingleFile, WdlValue}

import scala.annotation.tailrec
import scala.collection.mutable.ArrayBuffer
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Random, Try}

object LocalBackend {
  // Folders
  val CromwellExecutionDir = "cromwell-executions"

  // Files
  val ReturnCodeFile = "rc"
  val StdoutFile = "stdout"
  val StderrFile = "stderr"
  val ScriptFile = "script"

  // Other
  val DockerFlag = "docker"

  def props(task: TaskDescriptor): Props = Props(new LocalBackend(task))
}

/**
  * Executes a task in local computer through command line. It can be also executed in a Docker container.
  * @param task Task descriptor.
  */
class LocalBackend(task: TaskDescriptor) extends BackendActor with StrictLogging {

  import LocalBackend._

  implicit val timeout = Timeout(10 seconds)
  private val subscriptions = ArrayBuffer[Subscription[ActorRef]]()
  private var processAbortFunc: Option[() => Unit] = None

  val workingDir = task.workingDir
  val taskWorkingDir = task.name
  val shardId = Random.nextInt(Integer.MAX_VALUE).toString;
  val executionDir = Paths.get(CromwellExecutionDir, workingDir, taskWorkingDir, shardId).toAbsolutePath
  val stdout = Paths.get(executionDir.toString, StdoutFile)
  val stderr = Paths.get(executionDir.toString, StderrFile)
  val script = Paths.get(executionDir.toString, ScriptFile)
  val returnCode = Paths.get(executionDir.toString, ReturnCodeFile)
  val stdoutWriter = stdout.untailed
  val stderrTailed = stderr.tailed(100)
  val argv = Seq("/bin/bash", script.toString)
  val dockerImage = getDockerImage(task.runtimeAttributes)
  val expressionEval = new WorkflowEngineFunctions(executionDir)


  /**
    * Prepare the task and context for execution.
    */
  override def prepare(): Unit = {
    logger.debug(s"Creating execution folder: $executionDir")
    executionDir.toString.toFile.createIfNotExists(true)

    try {
      val command = initiateCommand()
      logger.debug(s"Creating bash script for executing command: ${command}.")
      writeBashScript(command, executionDir)
      notifyToSubscribers(new TaskStatus(Status.Created))
    } catch {
      case ex: Exception => notifyToSubscribers(new TaskFinalStatus(Status.Failed, FailureResult(ex)))
    }
  }

  /**
    * Stops a task execution.
    */
  override def stop(): Unit = {
    processAbortFunc.get.apply()
    notifyToSubscribers(new TaskStatus(Status.Canceled))
  }

  /**
    * Executes task in given context.
    */
  override def execute(): Unit = {
    notifyToSubscribers(executeTask)
  }

  /**
    * Performs a cleanUp after the task was executed.
    */
  override def cleanUp(): Unit = ()

  /**
    * Subscribe to events on backend.
    */
  override def subscribeToEvent[A](subscription: Subscription[A]): Unit = {
    val sub = subscription.asInstanceOf[Subscription[ActorRef]]
    subscriptions += sub
    sub.subscriber ! Subscribed
  }

  /**
    * Unsubscribe to events on backend.
    */
  override def unsubscribeToEvent[A](subscription: Subscription[A]): Unit = {
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
    * Notifies to subscribers about a new event while executing the task.
    * @param message A task status event.
    */
  private def notifyToSubscribers(message: ExecutionEvent): Unit = {
    subscriptions.filter(subs => subs.eventType.isInstanceOf[ExecutionEvent]).foreach(
      subs => subs.subscriber ! message)
  }

  /**
    * Gather Docker image name from runtime attributes. It it's not present returns none.
    * @param runtimeAttributes Runtime requirements for the task.
    */
  private def getDockerImage(runtimeAttributes: Map[String, String]): Option[String] = {
    if (runtimeAttributes.filter(p => p._1.equals(DockerFlag)).size == 1 &&
      !runtimeAttributes.filter(p => p._1.equals(DockerFlag)).values.head.isEmpty) {
      Some(runtimeAttributes.filter(p => p._1.equals(DockerFlag)).values.head)
    } else None
  }

  /**
    * Extracts folder path from specific input file.
    * @param file Absolute path from input file.
    * @return File's folder.
    */
  private def extractFolder(file: String): Option[String] = {
    try {
      Some(file.substring(0, file.lastIndexOf("/")))
    } catch {
      case oooe: StringIndexOutOfBoundsException =>
        logger.warn("Input with no valid folder pattern. It may be a intermediate value.", oooe)
        None
      case ex: Exception =>
        logger.error("Unhandled exception.", ex)
        throw ex
    }
  }

  /**
    * Creates docker command in order to execute the task into a container.
    * @param image Docker image name.
    * @return Command to execute.
    */
  private def buildDockerRunCommand(image: String): String = {
    val dockerVolume = "-v %s:%s"
    val inputFolder = task.inputs.filter(p => p._2.isInstanceOf[WdlFile]).map(p => extractFolder(p._2.toWdlString.replace("\"", "")))
      .filter(p => p.isDefined).map(p => p.get).toList.distinct
    val inputVolumes = inputFolder match {
      case a: List[String] => inputFolder.map(v => dockerVolume.format(v, v)).mkString(" ")
      case _ => ""
    }
    val outputVolume = dockerVolume.format(executionDir, executionDir)
    log.debug(s"DockerInputVolume: $inputVolumes")
    log.debug(s"DockerOutputVolume: $outputVolume")
    val dockerCmd = s"docker run %s %s --rm %s %s" //TODO: make it configurable from file.
    dockerCmd.format(inputVolumes, outputVolume, image, argv.mkString(" ")).replace("  ", " ")
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
    * Resolves absolute path for output files. If it is not a file, it will returns same value.
    * @param output Pair of WdlType and WdlValue
    * @return WdlValue with absolute path if it is a file.
    */
  private def resolveOutputValue(output: (WdlType, Try[WdlValue])): WdlValue = {
    def getAbsolutePath(file: String): Path = {
      val absolutePath = Paths.get(executionDir.toString, file)
      absolutePath.exists match {
        case true => absolutePath
        case false => throw new IllegalStateException(s"Output file $file does not exist.")
      }
    }

    //TODO: WdlArray[File] is missing here.
    output._1 match {
      case lhs if output._2.get.wdlType == WdlFileType | lhs == WdlFileType =>
        new WdlSingleFile(getAbsolutePath(output._2.get.toWdlString.replace("\"", "").trim).toString)
      case lhs if lhs == output._2.get.wdlType => output._2.get
      case lhs => lhs.coerceRawValue(output._2.get).get
    }
  }

  /**
    * Run a command using a bash script.
    * @return A TaskStatus with the final status of the task.
    */
  private def executeTask(): TaskFinalStatus = {
    def getCmdToExecute(): Seq[String] = dockerImage match {
      case Some(image) => buildDockerRunCommand(image).split(" ").toSeq
      case None => argv
    }

    val process = getCmdToExecute.run(ProcessLogger(stdoutWriter writeWithNewline, stderrTailed writeWithNewline))
    processAbortFunc = Some(() => process.destroy())
    notifyToSubscribers(new TaskStatus(Status.Running))
    val backendCommandString = argv.map(s => "\"" + s + "\"").mkString(" ")
    logger.debug(s"command: $backendCommandString")
    val processReturnCode = process.exitValue() // blocks until process finishes

    Vector(stdoutWriter.writer, stderrTailed.writer) foreach {
      _.flushAndClose()
    }

    val stderrFileLength = Try(Files.size(stderr)).getOrElse(0L)

    // TODO: check continue on error.

    if (stderrFileLength > 0) {
      TaskFinalStatus(Status.Failed, FailureTaskResult(
        new IllegalStateException("StdErr file is not empty."), processReturnCode, stderr.toString))
    } else if (stderrFileLength <= 0 && processReturnCode != 0) {
      TaskFinalStatus(Status.Failed, FailureTaskResult(
        new IllegalStateException("RC code is not equals to zero."), processReturnCode, stderr.toString))
    } else {
      def lookupFunction: String => WdlValue = WdlExpression.standardLookupFunction(task.inputs, task.declarations, expressionEval)
      val outputsExpressions = task.outputs.map(
        output => output.name ->(output.wdlType, output.expression.evaluate(lookupFunction, expressionEval)))
      processOutputResult(processReturnCode, outputsExpressions)
    }
  }

  /**
    * Process output evaluating expressions, checking for created files and converting WdlString to WdlSimpleFile if necessary.
    * @param processReturnCode Return code from process.
    * @param outputsExpressions Outputs.
    * @return TaskStatus with final status of the task.
    */
  private def processOutputResult(processReturnCode: Int, outputsExpressions: Seq[(String, (WdlType, Try[WdlValue]))]): TaskFinalStatus = {
    if (outputsExpressions.filter(_._2._2.isFailure).size > 0) {
      TaskFinalStatus(Status.Failed, FailureTaskResult(new IllegalStateException("Failed to evaluate output expressions.",
        outputsExpressions.filter(_._2._2.isFailure).head._2._2.failed.get), processReturnCode, stderr.toString))
    } else {
      try {
        TaskFinalStatus(Status.Succeeded, SuccessfulTaskResult(outputsExpressions.map(
          output => output._1 -> resolveOutputValue(output._2)).toMap))
      } catch {
        case ex: Exception => TaskFinalStatus(Status.Failed, FailureTaskResult(ex, processReturnCode, stderr.toString))
      }
    }
  }

  /**
    * 1) Remove all leading newline chars
    * 2) Remove all trailing newline AND whitespace chars
    * 3) Remove all *leading* whitespace that's common among every line in the input string
    * For example, the input string:
    * "
    * first line
    * second line
    * third line
    *
    * "
    * Would be normalized to:
    * "first line
    * second line
    * third line"
    * @param s String to process
    * @return String which has common leading whitespace removed from each line
    */
  private def normalize(s: String): String = {
    val trimmed = stripAll(s, "\r\n", "\r\n \t")
    val parts = trimmed.split("[\\r\\n]+")
    val indent = parts.map(leadingWhitespaceCount).min
    parts.map(_.substring(indent)).mkString("\n")
  }

  private def leadingWhitespaceCount(s: String): Int = {
    val Ws = Pattern.compile("[\\ \\t]+")
    val matcher = Ws.matcher(s)
    if (matcher.lookingAt) matcher.end else 0
  }

  private def stripAll(s: String, startChars: String, endChars: String): String = {
    /* https://stackoverflow.com/questions/17995260/trimming-strings-in-scala */
    @tailrec
    def start(n: Int): String = {
      if (n == s.length) ""
      else if (startChars.indexOf(s.charAt(n)) < 0) end(n, s.length)
      else start(1 + n)
    }

    @tailrec
    def end(a: Int, n: Int): String = {
      if (n <= a) s.substring(a, n)
      else if (endChars.indexOf(s.charAt(n - 1)) < 0) s.substring(a, n)
      else end(a, n - 1)
    }

    start(0)
  }

  private def initiateCommand(): String = {
    normalize(task.commandTemplate.map(_.instantiate(task.declarations, task.inputs, expressionEval)).mkString(""))
  }
}
