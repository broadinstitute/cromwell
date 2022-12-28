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
//import akka.util.Timeout
//import cromwell.backend.google.pipelines.common.api.RunStatus.TerminalRunStatus

object GcpBatchAsyncBackendJobExecutionActor {

  type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, GcpBatchRunStatus]

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams) extends BackendJobLifecycleActor with StandardAsyncExecutionActor with AskSupport with GcpBatchRunCreationClient with GcpBatchStatusRequestClient with CromwellInstrumentation {

  import GcpBatchAsyncBackendJobExecutionActor._

  lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id

  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = GcpBatchRunStatus

  // temporary until GCP Batch client can generate random job IDs
  val jobTemp = "job-" + java.util.UUID.randomUUID.toString

  //override def receive: Receive = pollingActorClientReceive orElse runCreationClientReceive orElse super.receive

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  override def dockerImageUsed: Option[String] = Option("test")

  //type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, StandardAsyncRunState]

  val backendSingletonActor: ActorRef = standardParams.backendSingletonActorOption.getOrElse(throw new RuntimeException("GCP Batch actor cannot exist without its backend singleton 2"))

  def uploadScriptFile(): Future[Unit] = {
    commandScriptContents
      .fold(
        errors => Future
          .failed(new RuntimeException(errors
            .toList
            .mkString(", "))),
        asyncIo
          .writeAsync(jobPaths
            .script, _, Seq
            .empty)
      )
  }
  // Primary entry point for cromwell to run GCP Batch job
  override def executeAsync(): Future[ExecutionHandle] = {

    val batchTest = BatchRequest(projectId="batch-testing-350715", region="us-central1", jobName=jobTemp)
    //val jobTest = GcpBatchJobName(jobName=jobTemp)

    for {
      _ <- uploadScriptFile()
      _ = backendSingletonActor ! batchTest
      //_ = backendSingletonActor ! jobTest
      runId = StandardAsyncJob(UUID.randomUUID().toString)  //temp to test


    } yield PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None)


   // val runId = StandardAsyncJob(UUID.randomUUID().toString)  //temp to test
    //log.info(s"runId ${runId}")
    //Future.successful(PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None))

  }

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(1
    .second, 5
    .minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 3
      .second, maxInterval = 20
      .second, multiplier = 1.1)

   override def pollStatusAsync(handle: GcpBatchPendingExecutionHandle): Future[StandardAsyncRunState] = {

     val jobId = handle.pendingJob.jobId
     val _ = handle
       .runInfo match {
       case Some(actualJob) =>
         actualJob
       case None =>
         throw new RuntimeException(
           s"pollStatusAsync called but job not available. This should not happen. Job Id $jobId"
         )
     }

     //implicit val timeout: Timeout = Timeout(5.seconds)

     for {

       answer <- super[GcpBatchStatusRequestClient].pollStatus(workflowId = workflowId, jobId = handle.pendingJob, gcpBatchJobId = jobTemp)

     }
       yield answer



     }

  override def isTerminal(runStatus: GcpBatchRunStatus): Boolean = {
    runStatus match {
      case _: GcpBatchRunStatus.TerminalRunStatus => true
      case _ => false
    }
  }

  override def isDone(runStatus: GcpBatchRunStatus): Boolean = {
    runStatus match {
      case _: GcpBatchRunStatus.Success =>
        println("GCP job matched isDone")
        true
      case _ => throw new RuntimeException(s"Cromwell programmer blunder: isSuccess was called on an incomplete RunStatus ($runStatus).")
    }
  }

  override val gcpBatchActor: ActorRef = backendSingletonActor

  override val gcpBatchApiActor: ActorRef = backendSingletonActor //added to test polling
}

