package cromwell.backend.impl.htcondor

import java.nio.file.{Files, Path}

import better.files._
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.impl.htcondor
import cromwell.core.PathFactory.{EnhancedPath, FlushingAndClosingWriter}
import cromwell.core.{TailedWriter, UntailedWriter}

import scala.annotation.tailrec
import scala.language.postfixOps
import scala.sys.process._

object JobStatus {
    val MapOfstatuses = Map(
      0 -> Created, // This is actually `unexpanded` in HtCondor, not sure what that actually means
      1 -> Created, // Idle
      2 -> Running,
      3 -> Removed,
      4 -> Completed,
      5 -> Failed, // SystemOnHold
      6 -> SubmissionError // Also the default
    )

  def fromCondorStatusCode(statusCode: Int): JobStatus = {
    MapOfstatuses.getOrElse(statusCode, SubmissionError) // By default we return SubmissionError
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

  def generateSubmitFile(path: Path, attributes: Map[String, String]): String = {
    def htCondorSubmitCommand(filePath: Path) = {
      s"${HtCondorCommands.Submit} ${filePath.toString}"
    }

    val submitFileWriter = path.untailed
    attributes.foreach(attribute => submitFileWriter.writeWithNewline(s"${attribute._1}=${attribute._2}"))
    submitFileWriter.writeWithNewline(HtCondorRuntimeKeys.Queue)
    submitFileWriter.writer.flushAndClose()
    logger.debug(s"submit file name is : $path")
    logger.debug(s"content of file is : ${path.lines.toList}")
    htCondorSubmitCommand(path)
  }

}

class HtCondorProcess extends StrictLogging {
  private val stdout = new StringBuilder
  private val stderr = new StringBuilder

  def processLogger: ProcessLogger = ProcessLogger(stdout append _, stderr append _)
  def processStdout: String = stdout.toString().trim
  def processStderr: String = stderr.toString().trim
  def commandList(command: String): Seq[String] = Seq("/bin/bash",command)
  def untailedWriter(path: Path): UntailedWriter = path.untailed
  def tailedWriter(limit: Int, path: Path): TailedWriter = path.tailed(limit)
  def externalProcess(cmdList: Seq[String], processLogger: ProcessLogger = processLogger): Process = cmdList.run(processLogger)

  /**
    * Returns the RC of this job when it finishes.  Sleeps and polls
    * until the 'rc' file is generated
    */
  def jobReturnCode(jobId: String, returnCodeFilePath: Path): Int = {

    @tailrec
    def recursiveWait(): Int =
      checkStatus(jobId) match {
        case status if JobStatus.isTerminal(status) =>
          Files.exists(returnCodeFilePath) match {
            case true => returnCodeFilePath.contentAsString.stripLineEnd.toInt
            case false =>
              val msg = s"JobStatus from Condor is terminal ($status) and no RC file exists!"
              logger.debug(msg)
              throw new IllegalStateException(msg)
          }
        case nonTerminalStatus =>
          Thread.sleep(5000)
          recursiveWait()
      }

    recursiveWait()
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
  val RequestMemory = "request_memory"
  val Cpu = "request_cpus"
  val Disk = "request_disk"
}
