package cromwell.backend.google.pipelines.batch

import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.backend._
import cromwell.core.WorkflowId

import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.async.ExecutionHandle
import akka.actor.ActorRef
import akka.pattern.AskSupport

import cromwell.services.instrumentation.CromwellInstrumentation

import java.util.UUID
import scala.concurrent.Future
import scala.concurrent.duration._

import GcpBatchBackendSingletonActor._

object GcpBatchAsyncBackendJobExecutionActor {

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams) extends BackendJobLifecycleActor with StandardAsyncExecutionActor with AskSupport with GcpBatchStatusRequestClient with CromwellInstrumentation {

  //import GcpBatchAsyncBackendJobExecutionActor._

  lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id

  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = GcpBatchRunStatus

  val jobTemp = "job-" + java.util.UUID.randomUUID.toString


  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */

  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  override def dockerImageUsed: Option[String] = Option("test")

  //type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, StandardAsyncRunState]

  val backendSingletonActor: ActorRef = standardParams.backendSingletonActorOption.getOrElse(throw new RuntimeException("GCP Batch actor cannot exist without its backend singleton 2"))

  // Primary entry point for cromwell to run GCP Batch job
  override def executeAsync(): Future[ExecutionHandle] = {

    val batchTest = BatchRequest(projectId="batch-testing-350715", region="us-central1", jobName=jobTemp)

    backendSingletonActor ! batchTest

    val runId = StandardAsyncJob(UUID.randomUUID().toString)  //temp to test
    log.info(s"runId ${runId}")
    Future.successful(PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None))

  }

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(1
    .second, 5
    .minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 3
      .second, maxInterval = 20
      .second, multiplier = 1.1)

   override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[StandardAsyncRunState] = {
     super[GcpBatchStatusRequestClient].pollStatus(workflowId = workflowId, jobId = handle.pendingJob, gcpBatchJobId = jobTemp)

     }

  override def isTerminal(runStatus: GcpBatchRunStatus): Boolean = ???

  override def isDone(runStatus: GcpBatchRunStatus): Boolean = {
    runStatus match {
      case _: GcpBatchRunStatus.Success => true
    }
  }

  override val gcpBatchActor: ActorRef = backendSingletonActor
}

