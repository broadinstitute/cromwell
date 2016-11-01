package cromwell.backend.impl.htcondor

import java.nio.file.{Files, Path}

import better.files._
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.impl.htcondor
import cromwell.core.path.{TailedWriter, UntailedWriter}
import cromwell.core.path.PathImplicits._
import cromwell.core.path.JavaWriterImplicits._

import scala.sys.process._

object JobStatus {
  val MapOfStatuses = Map(
    0 -> Created, // This is actually `unexpanded` in HtCondor, not sure what that actually means
    1 -> Created, // Idle
    2 -> Running,
    3 -> Removed,
    4 -> Completed,
    5 -> Failed, // SystemOnHold
    6 -> SubmissionError // Also the default
  )

  def fromCondorStatusCode(statusCode: Int): JobStatus = {
    MapOfStatuses.getOrElse(statusCode, SubmissionError) // By default we return SubmissionError
  }

  def isTerminal(jobStatus: JobStatus): Boolean = jobStatus.isInstanceOf[TerminalJobStatus]
}

sealed trait JobStatus
sealed trait TerminalJobStatus extends JobStatus
case object Created extends JobStatus
case object Running extends JobStatus
case object Completed extends TerminalJobStatus
case object Removed extends TerminalJobStatus
case object Failed extends TerminalJobStatus
case object Aborted extends TerminalJobStatus
case object SubmissionError extends TerminalJobStatus


object HtCondorCommands {
  val SubmitOutputPattern = "(\\d*) job\\(s\\) submitted to cluster (\\d*)\\."
  val Submit = "condor_submit"
  val Remove = "condor_rm"
  private val JobStatus = "condor_q %s -autoformat JobStatus"
  def generateJobStatusCommand(jobId: String): String = HtCondorCommands.JobStatus.format(jobId)
}

class HtCondorCommands extends StrictLogging {

  /**
    * Writes the script file containing the user's command from the WDL as well
    * as some extra shell code for monitoring jobs
    */
  def writeScript(instantiatedCommand: String, filePath: Path, containerRoot: Path): Unit = {
    logger.debug(s"Writing bash script for execution. Command: $instantiatedCommand.")
    val scriptBody = s"""

#!/bin/sh
cd $containerRoot
$instantiatedCommand
echo $$? > rc

""".trim + "\n"
    File(filePath).write(scriptBody)
    ()
  }

  def generateSubmitFile(path: Path, attributes: Map[String, Any], nativeSpecs: Option[Array[String]]): String = {
    def htCondorSubmitCommand(filePath: Path) = {
      s"${HtCondorCommands.Submit} ${filePath.toString}"
    }

    val submitFileWriter = path.untailed
    attributes.foreach { attribute => submitFileWriter.writeWithNewline(s"${attribute._1}=${attribute._2}") }
    //Native specs is intended for attaching HtCondor native configuration such as 'requirements' and 'rank' definition
    //directly to the submit file.
    nativeSpecs foreach { _.foreach { submitFileWriter.writeWithNewline } }
    submitFileWriter.writeWithNewline(HtCondorRuntimeKeys.Queue)
    submitFileWriter.writer.flushAndClose()
    logger.debug(s"submit file name is : $path")
    logger.debug(s"content of file is : ${File(path).lines.toList}")
    htCondorSubmitCommand(path)
  }

}

class HtCondorProcess extends StrictLogging {
  private val stdout = new StringBuilder
  private val stderr = new StringBuilder

  def processLogger: ProcessLogger = ProcessLogger(s => { stdout append s; () }, s => { stderr append s; () })
  def processStdout: String = stdout.toString().trim
  def processStderr: String = stderr.toString().trim
  def commandList(command: String): Seq[String] = Seq("/bin/bash",command)
  def untailedWriter(path: Path): UntailedWriter = path.untailed
  def tailedWriter(limit: Int, path: Path): TailedWriter = path.tailed(limit)
  def externalProcess(cmdList: Seq[String], processLogger: ProcessLogger = processLogger): Process = cmdList.run(processLogger)

  /**
    * Returns the RC of this job if it has finished.
    */
  def jobReturnCode(jobId: String, returnCodeFilePath: Path): Option[Int] = {

    checkStatus(jobId) match {
      case status if JobStatus.isTerminal(status) =>
        Files.exists(returnCodeFilePath) match {
          case true => Option(File(returnCodeFilePath).contentAsString.stripLineEnd.toInt)
          case false =>
            val msg = s"JobStatus from Condor is terminal ($status) and no RC file exists!"
            logger.debug(msg)
            throw new IllegalStateException(msg)
        }
      case nonTerminalStatus => None
    }
  }

  private def checkStatus(jobId: String): JobStatus = {
    val htCondorProcess = new HtCondorProcess
    val commandArgv = HtCondorCommands.generateJobStatusCommand(jobId).split(" ").toSeq
    val process = htCondorProcess.externalProcess(commandArgv)
    val returnCode = process.exitValue()
    returnCode match {
      case 0 =>
        val stdout = htCondorProcess.processStdout
        // If stdout is empty, that means the job got removed from the queue. Return Completed in that case
        val status = if (stdout.isEmpty) htcondor.Completed else JobStatus.fromCondorStatusCode(htCondorProcess.processStdout.toInt)
        logger.info("Condor JobId {} current status: {}", jobId, status)
        status
      case errorCode =>
        val msg = "Could not retreive status from the queue: " + htCondorProcess.processStderr
        logger.error(msg)
        throw new IllegalStateException(msg)
    }
  }

}

object HtCondorRuntimeKeys {
  val Executable = "executable"
  val Arguments = "arguments"
  val Error = "error"
  val Output = "output"
  val Log = "log"
  val Queue = "queue"
  val Rank = "rank"
  val Requirements = "requirements"
  val Cpu = "request_cpus"
  val Memory = "request_memory"
  val Disk = "request_disk"
  val LogXml = "log_xml"
  val LeaveInQueue = "leave_in_queue"
  val InitialWorkingDir = "Iwd"
}
