package cromwell.engine.backend.sge

import java.nio.file.Files

import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.{Call, CallInputs, CallOutputs, WorkflowDescriptor}
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.backend.local.{LocalBackend, SharedFileSystem}
import cromwell.engine.backend.{Backend, TaskAbortedException}
import cromwell.engine.db.DataAccess
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, _}
import cromwell.parser.BackendType
import cromwell.util.FileUtil._

import scala.annotation.tailrec
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

class SgeBackend extends Backend with SharedFileSystem with LazyLogging {
  type BackendCall = SgeBackendCall

  import LocalBackend.WriteWithNewline

  override def backendType = BackendType.SGE

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        key: CallKey,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortRegistrationFunction): BackendCall = {
    SgeBackendCall(this, workflowDescriptor, key, locallyQualifiedInputs, abortRegistrationFunction)
  }

  override def execute(backendCall: BackendCall): Try[CallOutputs] =  {
    val tag = makeTag(backendCall)
    backendCall.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"$tag `$instantiatedCommand`")
        writeScript(backendCall, instantiatedCommand)
        launchQsub(backendCall) match {
          case (qsubRc, _) if qsubRc != 0 => Failure(new Throwable(s"Error: qsub exited with return code: $qsubRc"))
          case (_, None) => Failure(new Throwable(s"Could not parse Job ID from qsub output"))
          case (_, Some(sgeJobId)) => pollForSgeJobCompletionThenPostProcess(backendCall, sgeJobId)
        }
      case Failure(ex) => Failure(ex)
    }
  }

  // TODO: Not much thought was given to this function
  override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)
                                 (implicit ec: ExecutionContext): Future[Any] = {
    Future.successful(Unit)
  }

  /**
   * Returns the RC of this job when it finishes.  Sleeps and polls
   * until the 'rc' file is generated
   */
  private def waitUntilComplete(backendCall: BackendCall): Int = {
    val tag = makeTag(backendCall)
    @tailrec
    def recursiveWait(): Int = Files.exists(backendCall.rc) match {
      case true => backendCall.rc.toFile.slurp.stripLineEnd.toInt
      case false =>
        logger.info(s"$tag 'rc' file does not exist yet")
        Thread.sleep(5000)
        recursiveWait()
    }
    recursiveWait()
  }

  /**
   * This returns a zero-parameter function which, when called, will kill the
   * SGE job with id `sgeJobId`.  It also writes to the 'rc' file the value '143'
   */
  private def killSgeJob(backendCall: BackendCall, sgeJobId: Int) = () => {
    val (_, qdelStdoutWriter) = backendCall.callRootPath.resolve("qdel.stdout").fileAndWriter
    val (_, qdelStderrWriter) = backendCall.callRootPath.resolve("qdel.stderr").fileAndWriter
    val argv = Seq("qdel", sgeJobId.toString)
    val process = argv.run(ProcessLogger(qdelStdoutWriter writeWithNewline, qdelStderrWriter writeWithNewline))
    val rc: Int = process.exitValue()
    Vector(qdelStdoutWriter, qdelStderrWriter) foreach { _.flushAndClose() }
    val (_, rcWriter) = backendCall.rc.fileAndWriter
    rcWriter.writeWithNewline("143")
    rcWriter.flushAndClose()
    logger.debug(s"qdel $sgeJobId (rc=$rc)")
  }

  /**
   * Writes the script file containing the user's command from the WDL as well
   * as some extra shell code for monitoring jobs
   */
  private def writeScript(backendCall: BackendCall, instantiatedCommand: String) = {
    val (_, scriptWriter) = backendCall.script.fileAndWriter
    scriptWriter.writeWithNewline("#!/bin/sh")
    scriptWriter.writeWithNewline(instantiatedCommand)
    scriptWriter.writeWithNewline("echo $? > rc")
    scriptWriter.flushAndClose()
  }

  /**
   * Launches the qsub command, returns a tuple: (rc, Option(sge_job_id))
   */
  private def launchQsub(backendCall: BackendCall): (Int, Option[Int]) = {
    val tag = makeTag(backendCall)
    val sgeJobName = s"cromwell_${backendCall.workflowDescriptor.shortId}_${backendCall.call.name}"
    val argv = Seq("qsub", "-terse", "-N", sgeJobName, "-V", "-b", "n", "-wd", backendCall.callRootPath.toAbsolutePath, "-o", backendCall.stdout.getFileName, "-e", backendCall.stderr.getFileName, backendCall.script.toAbsolutePath).map(_.toString)
    val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
    logger.info(s"$tag backend command: $backendCommandString")

    val (qsubStdoutFile, qsubStdoutWriter) = backendCall.callRootPath.resolve("qsub.stdout").fileAndWriter
    val (qsubStderrFile, qsubStderrWriter) = backendCall.callRootPath.resolve("qsub.stderr").fileAndWriter
    val process = argv.run(ProcessLogger(qsubStdoutWriter writeWithNewline, qsubStderrWriter writeWithNewline))
    val rc: Int = process.exitValue()
    Vector(qsubStdoutWriter, qsubStderrWriter) foreach { _.flushAndClose() }

    // The -terse option to qsub makes it so stdout only has the job ID, if it was successfully launched
    val jobId = Try(qsubStdoutFile.slurp.stripLineEnd.toInt) match {
      case Success(id) => Some(id)
      case Failure(ex) =>
        logger.error(s"$tag Could not find SGE job ID from qsub stdout file.\n\nCheck the qsub stderr file for possible errors: ${qsubStderrFile.toAbsolutePath.toString}")
        None
    }
    logger.info(s"$tag qsub rc=$rc, job ID=${jobId.getOrElse("NONE")}")
    (rc, jobId)
  }

  /**
   * This waits for a given SGE job to finish.  When finished, it post-processes the job
   * and returns the outputs for the call
   */
  private def pollForSgeJobCompletionThenPostProcess(backendCall: BackendCall, sgeJobId: Int): Try[CallOutputs] = {
    val tag = makeTag(backendCall)
    val abortFunction = killSgeJob(backendCall, sgeJobId)
    val waitUntilCompleteFunction = waitUntilComplete(backendCall)
    backendCall.callAbortRegistrationFunction.register(AbortFunction(abortFunction))
    val jobRc = waitUntilCompleteFunction
    logger.info(s"$tag SGE job completed (rc=$jobRc)")
    (jobRc, backendCall.stderr.toFile.length) match {
      case (r, _) if r == 143 =>
        Failure(new TaskAbortedException()) // Special case to check for SIGTERM exit code - implying abort
      case (r, _) if r != 0 =>
        Failure(new Throwable(s"$tag SGE job failed because of non-zero return code: $r"))
      case (_, stderrLength) if stderrLength > 0 && backendCall.call.failOnStderr =>
        Failure(new Throwable(s"$tag SGE job failed because there were $stderrLength bytes on standard error"))
      case (_, _) =>
        postProcess(backendCall)
    }
  }
}
