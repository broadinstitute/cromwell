package cromwell.engine.backend.sge

import java.io.File
import java.nio.file.{Files, Path, Paths}

import com.typesafe.scalalogging.LazyLogging
import cromwell.binding.types.{WdlFileType, WdlMapType}
import cromwell.binding.values.{WdlValue, _}
import cromwell.binding.{Call, CallInputs, CallOutputs, WorkflowDescriptor, _}
import cromwell.engine.backend.Backend.RestartableWorkflow
import cromwell.engine.db.{CallStatus, DataAccess}
import cromwell.engine.{AbortFunctionRegistration, _}
import cromwell.engine.backend.local.{LocalFileSystemOperations, LocalBackend}
import cromwell.engine.backend.{Backend, BackendCall, TaskAbortedException}
import cromwell.parser.BackendType
import cromwell.util.FileUtil._
import scala.concurrent.{ExecutionContext, Future}

import scala.annotation.tailrec
import scala.language.postfixOps
import scala.sys.process._
import scala.util.{Failure, Success, Try}

class SgeBackend extends Backend with LocalFileSystemOperations with LazyLogging {
  type T = SgeBackendCall

  import LocalBackend.WriteWithNewline

  override def bindCall(workflowDescriptor: WorkflowDescriptor,
                        call: Call,
                        locallyQualifiedInputs: CallInputs,
                        abortRegistrationFunction: AbortFunctionRegistration): BackendCall = {
    SgeBackendCall(this, workflowDescriptor, call, locallyQualifiedInputs, abortRegistrationFunction)
  }

  override def execute(backendCall: T): Try[CallOutputs] =  {
    val tag: String = s"${this.getClass.getName} [UUID(${backendCall.workflowDescriptor.shortId}):${backendCall.call.name}]"

    backendCall.instantiateCommand match {
      case Success(instantiatedCommand) =>
        logger.info(s"$tag `$instantiatedCommand`")

        val (scriptFile, scriptWriter) = backendCall.script.fileAndWriter
        scriptWriter.writeWithNewline("#!/bin/sh")
        scriptWriter.writeWithNewline(instantiatedCommand)
        scriptWriter.writeWithNewline("echo $? > rc")
        scriptWriter.flushAndClose()

        val sgeJobName = s"cromwell_${backendCall.workflowDescriptor.shortId}_${backendCall.call.name}"
        val argv = Seq("qsub", "-N", sgeJobName, "-V", "-b", "n", "-wd", backendCall.callRootPath.toAbsolutePath, "-o", backendCall.stdout.getFileName, "-e", backendCall.stderr.getFileName, backendCall.script.toAbsolutePath).map(_.toString)
        val backendCommandString = argv.map(s => "\""+s+"\"").mkString(" ")
        logger.info(s"$tag backend command: $backendCommandString")

        val (qsubStdoutFile, qsubStdoutWriter) = backendCall.callRootPath.resolve("qsub.stdout").fileAndWriter
        val (qsubStderrFile, qsubStderrWriter) = backendCall.callRootPath.resolve("qsub.stderr").fileAndWriter
        val process = argv.run(ProcessLogger(qsubStdoutWriter writeWithNewline, qsubStderrWriter writeWithNewline))
        val rc: Int = process.exitValue()
        Vector(qsubStdoutWriter, qsubStderrWriter) foreach { _.flushAndClose() }
        logger.info(s"$tag qsub rc=$rc")

        "^Your job (\\d+)".r.findFirstMatchIn(qsubStdoutFile.slurp) match {
          case Some(m) =>
            val sgeJobId = m.group(1).toInt

            def killSgeJob(): Unit = {
              val (qdelStdoutFile, qdelStdoutWriter) = backendCall.callRootPath.resolve("qdel.stdout").fileAndWriter
              val (qdelStderrFile, qdelStderrWriter) = backendCall.callRootPath.resolve("qdel.stderr").fileAndWriter
              val argv = Seq("qdel", sgeJobId.toString)
              val process = argv.run(ProcessLogger(qdelStdoutWriter writeWithNewline, qdelStderrWriter writeWithNewline))
              val rc: Int = process.exitValue()
              Vector(qdelStdoutWriter, qdelStderrWriter) foreach { _.flushAndClose() }
              logger.debug(s"qsub $sgeJobId (rc=$rc)")
            }

            backendCall.callAbortRegistrationFunction.register(AbortFunction(killSgeJob))

            @tailrec
            def waitUntilComplete(): Int = Files.exists(backendCall.rc) match {
              case true => backendCall.rc.toFile.slurp.stripLineEnd.toInt
              case false =>
                logger.info(s"$tag File ${backendCall.rc} does not exist yet")
                Thread.sleep(5000)
                waitUntilComplete()
            }

            val jobRc = waitUntilComplete()

            logger.info(s"$tag File ${backendCall.rc} now exists")

            val stderrFileLength = backendCall.stderr.toFile.length

            if (backendCall.call.failOnStderr && stderrFileLength > 0) {
              logger.info(s"$tag stderr fail")
              Failure(new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: stderr has length $stderrFileLength for command: $instantiatedCommand"))
            } else {
              logger.info(s"$tag jobRc=$jobRc")
              jobRc match {
                case 0 =>
                  logger.info(s"$tag postProcess")
                  postProcess(backendCall)
                case 143 => Failure(new TaskAbortedException()) // Special case to check for SIGTERM exit code - implying abort
                case _ => Failure(new Throwable(s"Workflow ${backendCall.workflowDescriptor.id}: return code $rc for command: $instantiatedCommand\n\nFull command was: ${argv.mkString(" ")}"))
              }
            }

          case None => Failure(new Throwable("Failed to find SGE job ID in qsub.stdout"))
        }
      case Failure(ex) => Failure(ex)
    }
  }

  // TODO: Not much thought was given to this function
  override def handleCallRestarts(restartableWorkflows: Seq[RestartableWorkflow], dataAccess: DataAccess)
                                 (implicit ec: ExecutionContext): Future[Any] = {
    Future.successful(Unit)
  }

  override def backendType = BackendType.SGE
}
