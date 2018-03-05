/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

package cromwell.backend.impl.aws

import software.amazon.awssdk.services.batch.model.JobStatus
import cromwell.core.ExecutionEvent
import scala.util.{Failure, Success, Try}
sealed trait RunStatus {
  import RunStatus._

  // Could be defined as false for Initializing and true otherwise, but this is more defensive.
  def isRunningOrComplete = this match {
    case Running | _: TerminalRunStatus => true
    case _ => false
  }
}

object RunStatus {
  def fromJobStatus(status: JobStatus, jobId: String, errorMessage: Option[String] = None,
                    eventList: Seq[ExecutionEvent] = Seq.empty): Try[RunStatus] = {
    status match {
      case JobStatus.FAILED => Success(Failed(jobId, errorMessage, eventList))
      case JobStatus.PENDING => Success(Initializing)
      case JobStatus.RUNNABLE => Success(Initializing)
      case JobStatus.RUNNING => Success(Running)
      case JobStatus.STARTING => Success(Running)
      case JobStatus.SUBMITTED => Success(Initializing)
      case JobStatus.SUCCEEDED => Success(Succeeded(eventList))
      // JobStatus.UNKNOWN_TO_SDK_VERSION
      case _ => Failure(new RuntimeException(s"job {$jobId} has an unknown status {$status}"))
    }
  }

  case object Initializing extends RunStatus
  case object Running extends RunStatus

  sealed trait TerminalRunStatus extends RunStatus {
    def eventList: Seq[ExecutionEvent]
  }

  sealed trait UnsuccessfulRunStatus extends TerminalRunStatus {
    val errorMessage: Option[String]
    lazy val prettyPrintedError: String = errorMessage map { e => s" Message: $e" } getOrElse ""
  }

  case class Succeeded(eventList: Seq[ExecutionEvent]) extends TerminalRunStatus {
    override def toString = "Succeeded"
  }

  object UnsuccessfulRunStatus {
    def apply(jobId: String, status: String, errorMessage: Option[String], eventList: Seq[ExecutionEvent]): UnsuccessfulRunStatus = {
      if (status == "Stopped") { // TODO: Verify this
        Stopped(jobId, errorMessage, eventList)
      } else {
        Failed(jobId, errorMessage, eventList)
      }
    }
  }

  final case class Stopped(jobId: String,
                          errorMessage: Option[String],
                          eventList: Seq[ExecutionEvent],
                          ) extends UnsuccessfulRunStatus {
    override def toString = "Stopped"
  }

  final case class Failed(jobId: String,
                          errorMessage: Option[String],
                          eventList: Seq[ExecutionEvent],
                          ) extends UnsuccessfulRunStatus {
    override def toString = "Failed"
  }

  final case class Cancelled(jobId: String,
                          errorMessage: Option[String],
                          eventList: Seq[ExecutionEvent],
                          ) extends UnsuccessfulRunStatus {
    override def toString = "Cancelled"
  }
}
