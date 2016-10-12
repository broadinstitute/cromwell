package cromwell.backend.sfs

import java.nio.file.Path

import akka.actor.{ActorRef, LoggingFSM}
import better.files.File
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.JobId
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.impl.sfs.config.ConfigInitializationData
import cromwell.backend.io.JobPaths
import cromwell.backend.validation.ValidatedRuntimeAttributes
import wdl4s.Call

import scala.concurrent.Promise
import scala.concurrent.duration._
import scala.util.matching.Regex

class ConfigBatchSubmissionActor(override val params: SharedFileSystemBatchSubmissionActorParams)
  extends SharedFileSystemBatchSubmissionActor {

  lazy val configInitializationData: ConfigInitializationData = params.backendInitializationDataOption match {
    case Some(data: ConfigInitializationData) => data
    case other => throw new RuntimeException(s"Unable to get config initialization data from $other")
  }

  override def getJob(exitValue: Int, stdout: Path, stderr: Path): SharedFileSystemJob = {
    val arrayJobIdRegex = configurationDescriptor.backendConfig.getString(ArrayJobIdRegexConfig).r
    val output = File(stdout).contentAsString.stripLineEnd
    output match {
      case arrayJobIdRegex(jobId) => SharedFileSystemJob(jobId)
      case _ =>
        throw new RuntimeException("Could not find job ID from stdout file. " +
          s"Check the stderr file for possible errors: $stderr")
    }
  }

  // TODO: Generic mapping of job to job+index
  override def jobIdWithIndex(jobId: JobId, index: Int) = ???
}

case class SharedFileSystemBatchSubmissionActorParams
(
  configurationDescriptor: BackendConfigurationDescriptor,
  backendInitializationDataOption: Option[BackendInitializationData]
)

trait SharedFileSystemBatchSubmissionActor extends BatchSubmissionActor {
  type EntryWithIndex = (BatchSubmissionEntry, Int)

  val params: SharedFileSystemBatchSubmissionActorParams

  lazy val configurationDescriptor = params.configurationDescriptor

  override def batchSubmitJobs(batchSubmissionData: BatchSubmissionData): Unit = {
    val grouped = batchSubmissionData.submissionQueue groupBy { entry =>
      SharedFileSystemBatchJobKey(entry.jobDescriptor.call, batchJobData(entry.jobData).validatedRuntimeAttributes)
    }

    grouped foreach {
      case (sharedFileSystemBatchJobKey, batchSubmissionEntries) =>
        submitBatch(sharedFileSystemBatchJobKey, batchSubmissionEntries, batchSubmissionData.aborting)
    }
  }

  // TODO: Make an implicit class?
  def batchJobData(jobData: BackendJobData): SharedFileSystemBatchJobData = {
    jobData.asInstanceOf[SharedFileSystemBatchJobData]
  }

  private def submitBatch(sharedFileSystemBatchJobKey: SharedFileSystemBatchJobKey,
                          batchSubmissionEntries: Seq[BatchSubmissionEntry], aborting: Seq[ActorRef]): Unit = {

    if (aborting.exists(abortingSender => batchSubmissionEntries.exists(_.sender == abortingSender))) {
      // Reply to all senders that the entry was aborted
      batchSubmissionEntries foreach (_.sender ! BatchSubmissionAborted)
    } else {
      submitBatch(sharedFileSystemBatchJobKey, batchSubmissionEntries)
    }
  }

  private def submitBatch(sharedFileSystemBatchJobKey: SharedFileSystemBatchJobKey,
                          batchSubmissionEntries: Seq[BatchSubmissionEntry]): Unit = {

    val entriesWithIndex: Seq[(BatchSubmissionEntry, Int)] = batchSubmissionEntries.zipWithIndex

    val taskScripts = entriesWithIndex map innerScript

    val taskIdVar = ArrayTaskVariableConfig

    val outerScript =
      s"""#!/bin/bash
          |
          |case $taskIdVar in
          |${taskScripts.mkString("\n")}
          |* )
          |  echo "Task index '$taskIdVar' is not recognized."
          |  exit 1
          |  ;;
          |esac
          |""".stripMargin

    def innerScript(entryWithIndex: EntryWithIndex): String = {
      val (batchSubmissionEntry, index) = entryWithIndex
      val taskIndex = index + 1
      val jobPaths = batchJobData(batchSubmissionEntry.jobData).jobPaths
      s"""$taskIndex )
         |  cd '${jobPaths.callExecutionRoot}' \\
         |  && /bin/bash '${jobPaths.script}' \\
         |  > '${jobPaths.stdout}' \\
         |  2> '${jobPaths.stderr}'
         |  ;;
         |""".stripMargin
    }

    val runner = new ProcessRunner(s"qsub -t 1-$count -V -b n $outerScriptPath $qsubVariablesFromValidatedRuntimeAttributes", "batch.submit", "batch.stderr")
    val exitValue = runner.run()
    if (exitValue != 0) {
      val failure = BatchSubmissionFailure(new RuntimeException("Unable to start job. " +
        s"Check the stderr file for possible errors: ${runner.stderrPath}"))
      batchSubmissionEntries foreach (_.sender ! failure)
    } else {
      val runningJob = getJob(exitValue, runner.stdoutPath, runner.stderrPath)
      tellJobIds(runningJob, entriesWithIndex)
    }

  }

  def getJob(exitValue: Int, stdout: Path, stderr: Path): SharedFileSystemJob

  // TODO: Generic mapping of job to job+index
  def jobIdWithIndex(jobId: JobId, index: Int): JobId

  def tellJobIds(jobId: JobId, entriesWithIndex: Seq[(BatchSubmissionEntry, Int)]): Unit = {
    entriesWithIndex foreach tellJobId(jobId)
  }

  def tellJobId(jobId: JobId)(entryWithIndex: EntryWithIndex): Unit = {
    val (batchSubmissionEntry, index) = entryWithIndex
    batchSubmissionEntry.sender ! BatchSubmissionSuccess(jobIdWithIndex(jobId, index))
  }
}

case class SharedFileSystemBatchJobData(validatedRuntimeAttributes: ValidatedRuntimeAttributes,
                                        jobPaths: JobPaths) extends BackendJobData

case class SharedFileSystemBatchJobKey(call: Call, validatedRuntimeAttributes: ValidatedRuntimeAttributes)

abstract class BatchSubmissionActor extends LoggingFSM[BatchSubmissionState, BatchSubmissionData] {

  startWith(BatchSubmissionIdle, BatchSubmissionData(Seq.empty, Seq.empty))

  when(BatchSubmissionQueueing, stateTimeout = 1.second) {
    case (StateTimeout, batchSubmissionData) => goto(BatchSubmissionIdle) using
      BatchSubmissionData(Seq.empty, Seq.empty)
  }

  whenUnhandled {
    // common code for both states
    case Event(entry: BatchSubmissionEntry, batchSubmissionData: BatchSubmissionData) =>
      goto(BatchSubmissionQueueing) using
        batchSubmissionData.copy(submissionQueue = batchSubmissionData.submissionQueue :+ entry)

    case Event(AbortJobCommand, batchSubmissionData: BatchSubmissionData) =>
      stay using batchSubmissionData.copy(aborting = batchSubmissionData.aborting :+ sender())

    case Event(fsmEvent, fsmState) =>
      log.warning("received unhandled request {} in state {}/{}", fsmEvent, stateName, fsmState)
      stay
  }

  onTransition {
    case BatchSubmissionQueueing -> BatchSubmissionIdle => batchSubmitJobs(stateData)
  }

  abstract def batchSubmitJobs(batchSubmissionData: BatchSubmissionData): Unit
}

sealed trait BatchSubmissionState

case object BatchSubmissionIdle extends BatchSubmissionState

case object BatchSubmissionQueueing extends BatchSubmissionState

trait BackendJobData

case class BatchSubmissionEntry(sender: ActorRef, jobDescriptor: BackendJobDescriptor, jobData: BackendJobData)

case class BatchSubmissionData(submissionQueue: Seq[BatchSubmissionEntry], aborting: Seq[ActorRef])

case object BatchSubmissionAborted

case class BatchSubmissionSuccess(jobId: JobId)

case class BatchSubmissionFailure(failure: Exception)
