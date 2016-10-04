package cromwell.backend.sfs

import akka.actor.{ActorRef, LoggingFSM}
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.validation.ValidatedRuntimeAttributes
import wdl4s.Call

import scala.concurrent.duration._


class SharedFileSystemBatchSubmissionActor extends BatchSubmissionActor {
  override def batchSubmitJobs(batchSubmissionData: BatchSubmissionData): Unit = {
    // TODO: Submit jobs and reply with job ids.

    val keyed = batchSubmissionData.queue groupBy { entry =>
      entry.jobData match {
        case SharedFileSystemBatchJobData(validatedRuntimeAttributes) =>
          SharedFileSystemBatchJobKey(entry.jobDescriptor.call, validatedRuntimeAttributes)
        case unexpected => throw new RuntimeException(s"Expected a SharedFileSystemBatchJobData. Got $unexpected")
      }
    }

    val outerScript =
      s"""
        |#!/bin/bash
        |
        |
        |
        |/bin/bash scripts/$$TASK_ID
      """.stripMargin

    val innerScript =
      s"""
        |$taskId) cd '$cwd' && /bin/bash '$script' > '$stdout' 2> '$stderr' ;;
      """.stripMargin

  }
}

case class SharedFileSystemBatchJobData(validatedRuntimeAttributes: ValidatedRuntimeAttributes) extends BackendJobData

case class SharedFileSystemBatchJobKey(call: Call, validatedRuntimeAttributes: ValidatedRuntimeAttributes)

abstract class BatchSubmissionActor extends LoggingFSM[BatchSubmissionState, BatchSubmissionData] {

  startWith(BatchSubmissionIdle, BatchSubmissionData(Seq.empty))

  when(BatchSubmissionQueueing, stateTimeout = 1.second) {
    case (StateTimeout, batchSubmissionData) => goto(BatchSubmissionIdle)
  }

  whenUnhandled {
    // common code for both states
    case Event(entry: BatchSubmissionEntry, batchSubmissionData: BatchSubmissionData) =>
      goto(BatchSubmissionQueueing) using batchSubmissionData.copy(queue = batchSubmissionData.queue :+ entry)

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

case class BatchSubmissionData(queue: Seq[BatchSubmissionEntry])
