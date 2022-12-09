package cromwell.backend.google.pipelines.batch

import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration._
import cromwell.backend._
import cromwell.backend.google.pipelines.common.Run
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.async.ExecutionHandle
import akka.actor.ActorRef

import java.time.OffsetDateTime
//import com.google.cloud.batch.v1.{BatchServiceClient, JobName}
import cromwell.services.instrumentation.CromwellInstrumentation

import java.util.UUID
import scala.concurrent.Future

object GcpBatchAsyncBackendJobExecutionActor {

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams) extends BackendJobLifecycleActor with StandardAsyncExecutionActor with CromwellInstrumentation {

  //import GcpBatchAsyncBackendJobExecutionActor._

  //override lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id

  /** The type of the run info when a job is started. */
  //override type StandardAsyncRunInfo = this.type
  override type StandardAsyncRunInfo = String
  //override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  //override type StandardAsyncRunState = this.type
  override type StandardAsyncRunState = String
  //override type StandardAsyncRunState = RunStatus


  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  //override def statusEquivalentTo(thiz: GcpBatchAsyncBackendJobExecutionActor
  //  .this.type)(that: GcpBatchAsyncBackendJobExecutionActor
  //  .this.type): Boolean = true

  //override def isTerminal(runStatus: GcpBatchAsyncBackendJobExecutionActor
  //  .this.type): Boolean = true

  override def statusEquivalentTo(thiz: String)(that: String): Boolean = thiz == that

  override def dockerImageUsed: Option[String] = Option("test")
  override def isTerminal(runStatus: String): Boolean = runStatus == "DummyDone"

  //type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, StandardAsyncRunState]

  val backendSingletonActor: ActorRef = standardParams.backendSingletonActorOption.getOrElse(throw new RuntimeException("GCP Batch actor cannot exist without its backend singleton 2"))

  // Primary entry point for cromwell to run GCP Batch job
  override def executeAsync(): Future[ExecutionHandle] = {

      val batchTest = BatchRequest(projectId="batch-testing-350715", region="us-central1", jobName="executor-test-1")

      backendSingletonActor ! batchTest
      val runId = StandardAsyncJob(UUID.randomUUID().toString)  //temp to test

      // jobDescriptor: BackendJobDescriptor,
      // pendingJob: BackendJobId,
      // runInfo: Option[BackendRunInfo],
      // previousState: Option[BackendRunState]

      Future.successful(PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None))
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

   //override def pollStatusAsync(handle: GcpBatchPendingExecutionHandle): Future[StandardAsyncRunState] = {
   //def pollStatusAsync(handle: GcpBatchPendingExecutionHandle): Future[String] = {
   override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[String] = {
     finishTime match {
       case Some(_) => Future.successful("running")
       case None => Future.successful("fake running")
       //case None => throw new RuntimeException("Status Async called but job not available job Id")
     }
    //val jobId = handle.pendingJob.jobId
    //val job = handle.runInfo match {
    //  case Some(actualJob) => actualJob
    //  case None => throw new RuntimeException(s"Status Async called but job not available job Id $jobId")
    //}
     // Future["test"]

    //val projectId = "batch-testing-350715"
    //val region = "us-central1"
    //val jobName = "test3"

    //val batchServiceClient = BatchServiceClient
    //  .create()

    //val job2 = batchServiceClient
    //  .getJob(JobName
    //    .newBuilder()
    //    .setProject(projectId)
    //    .setLocation(region)
    //    .setJob(jobName)
    //    .build())

    //job2
  }
}

case class BatchRequest(projectId: String, region: String, jobName: String)
