package cromwell.backend.google.pipelines.batch

import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.backend._
import cromwell.core.WorkflowId

//import scala.concurrent.Promise
//import cromwell.backend.google.pipelines.common.api.RunStatus
//import cromwell.backend.google.pipelines.common.Run
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.async.ExecutionHandle
import akka.actor.ActorRef
import akka.pattern.AskSupport

import java.time.OffsetDateTime
import cromwell.services.instrumentation.CromwellInstrumentation

import java.util.UUID
import scala.concurrent.Future
import scala.concurrent.duration._

object GcpBatchAsyncBackendJobExecutionActor {

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams) extends BackendJobLifecycleActor with StandardAsyncExecutionActor with AskSupport with GcpBatchStatusRequestClient with CromwellInstrumentation {

  //import GcpBatchAsyncBackendJobExecutionActor._



  lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id

  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = GcpBatchRunStatus


  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */

  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  override def dockerImageUsed: Option[String] = Option("test")

  //override def isTerminal(runStatus: String): Boolean = runStatus == "DummyDone"
  override def isTerminal(runStatus: GcpBatchRunStatus): Boolean = ???

  //type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, StandardAsyncRunState]

  val backendSingletonActor: ActorRef = standardParams.backendSingletonActorOption.getOrElse(throw new RuntimeException("GCP Batch actor cannot exist without its backend singleton 2"))

  // Primary entry point for cromwell to run GCP Batch job
  override def executeAsync(): Future[ExecutionHandle] = {

    val batchTest = BatchRequest(projectId="batch-testing-350715", region="us-central1", jobName="executor-test-2")

    backendSingletonActor ! batchTest

    val runId = StandardAsyncJob(UUID.randomUUID().toString)  //temp to test
    Future.successful(PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None))

    //val runBatchResponse = for {
      //completionPromise = Promise
      //_ = backendSingletonActor ! batchTest

    //} yield batchRunId

     // runBatchResponse map { batchRunId =>PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None) }
      //}



      // jobDescriptor: BackendJobDescriptor,
      // pendingJob: BackendJobId,
      // runInfo: Option[BackendRunInfo],
      // previousState: Option[BackendRunState]

     // Future.successful(PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None))
      //Future.successful(
      //  PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](
      //    jobDescriptor = jobDescriptor,
      //    pendingJob = StandardAsyncJob(UUID.randomUUID().toString),
      //    runInfo = Option[StandardAsyncRunInfo],
      //    previousState = None
      //  )
      //)
  }

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(1
    .second, 5
    .minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 3
      .second, maxInterval = 20
      .second, multiplier = 1.1)

  var finishTime: Option[OffsetDateTime] = None


   //def pollStatusAsync(handle: GcpBatchPendingExecutionHandle): Future[String] = {
   override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[StandardAsyncRunState] = {
     super[GcpBatchStatusRequestClient].pollStatus(workflowId = workflowId, jobId = handle.pendingJob)

     }

  override val gcpBatchActor: ActorRef = backendSingletonActor
}

case class BatchRequest(projectId: String, region: String, jobName: String)
