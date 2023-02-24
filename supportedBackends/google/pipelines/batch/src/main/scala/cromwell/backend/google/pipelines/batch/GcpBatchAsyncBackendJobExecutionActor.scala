package cromwell.backend.google.pipelines.batch

import akka.util.Timeout
import com.google.api.gax.rpc.NotFoundException
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.backend._
import java.util.concurrent.ExecutionException
import scala.concurrent.Await
import cromwell.core.{ExecutionEvent, WorkflowId}
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.async.ExecutionHandle
import akka.actor.ActorRef
import akka.pattern.AskSupport
import cromwell.services.instrumentation.CromwellInstrumentation
import scala.concurrent.Future
import scala.concurrent.duration._
import GcpBatchBackendSingletonActor._
import cromwell.backend.google.pipelines.batch.RunStatus.{Running, Succeeded, TerminalRunStatus}

object GcpBatchAsyncBackendJobExecutionActor {

  type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, RunStatus]

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor
    with StandardAsyncExecutionActor
    with AskSupport
    with GcpBatchJobCachingActorHelper
    with GcpBatchStatusRequestClient
    with CromwellInstrumentation {

  import GcpBatchAsyncBackendJobExecutionActor._

  lazy val gcpBatchCommand: String = jobDescriptor.taskCall.callable.commandTemplateString(Map.empty)

  lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id

  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = RunStatus

  // temporary until GCP Batch client can generate random job IDs
  private val jobTemp = "job-" + java.util.UUID.randomUUID.toString

  //override def receive: Receive = pollingActorClientReceive orElse runCreationClientReceive orElse super.receive
  //override def receive: Receive = pollingActorClientReceive orElse super.receive

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  private lazy val jobDockerImage = jobDescriptor.maybeCallCachingEligible.dockerHash
                                                 .getOrElse(runtimeAttributes.dockerImage)

  override def dockerImageUsed: Option[String] = Option(jobDockerImage)

  //type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, StandardAsyncRunState]

  val backendSingletonActor: ActorRef = standardParams.backendSingletonActorOption
                                                      .getOrElse(throw new RuntimeException("GCP Batch actor cannot exist without its backend singleton 2"))

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

    val runBatchResponse = for {
      _ <- uploadScriptFile()
      _ = backendSingletonActor ! BatchRequest(workflowId, projectId = "batch-testing-350715", region = "us-central1", jobName = jobTemp, runtimeAttributes, gcpBatchCommand)
      runId = StandardAsyncJob(jobTemp)

    }
    yield runId

    runBatchResponse map { runId => PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None) }

  }

  override def reconnectAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = {
    log.info("reconnect async runs") // in for debugging remove later
    val handle = PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](jobDescriptor, jobId, Option(Run(jobId)), previousState = None)
    Future.successful(handle)
  }

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(5
    .second, 5
    .minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 5
      .seconds, maxInterval = 20
      .seconds, multiplier = 1.1)

  override def pollStatusAsync(handle: GcpBatchPendingExecutionHandle): Future[RunStatus] = {

    val jobId = handle.pendingJob.jobId

    val job = handle.runInfo match {
      case Some(actualJob) => actualJob
      case None =>
        throw new RuntimeException(
          s"pollStatusAsync called but job not available. This should not happen. Job Id $jobId"
        )
    }

    log.info(s"started polling for ${job} with jobId ${jobId}")

    //super[GcpBatchStatusRequestClient].pollStatus(workflowId, handle.pendingJob, jobTemp)

    //temporary. Added to resolve issue with poll async starter before job submitted.
    implicit val timeout: Timeout = Timeout(60.seconds)
    val futureResult = backendSingletonActor ? BatchJobAsk(jobId)
    val result = Await.result(futureResult, timeout.duration).asInstanceOf[String]
    log.info(result)

    try {
      val gcpBatchPoll = new GcpBatchJobGetRequest
      val result = gcpBatchPoll.GetJob(jobId)
      RunStatus.fromJobStatus(result)

    }
    catch {
      case nfe: NotFoundException => //added to account for job not found errors because polling async happens before job is submitted
        nfe.printStackTrace()
        Future.successful(Running)
      case ee: ExecutionException => //added to account for job not found errors because polling async happens before job is submitted
        ee.printStackTrace()
        Future.successful(Running)
    }

  }

  override def isTerminal(runStatus: RunStatus): Boolean = {
    //runStatus.isTerminal

    runStatus match {
      case jobSucceeded: RunStatus.Succeeded =>
        log.info("isTerminal match Succeeded running with status {}", jobSucceeded)
        true
      case jobFailed: RunStatus.Failed =>
        log.info("isTerminal match Failed with status {}", jobFailed)
        false
      case _: TerminalRunStatus =>
        val tempTermStatus = runStatus.toString
        log.info(f"isTerminal match TerminalRunStatus running with status $tempTermStatus")
        true
      case other =>
        log.info(f"isTerminal match _ running with status $other")
        false
    }
  }

  override def isDone(runStatus: RunStatus): Boolean = {
    runStatus match {
      case _: RunStatus.Succeeded =>
        log.info("GCP job matched isDone")
        true
      case _ =>
        log.info("did not match isDone")
        false //throw new RuntimeException(s"Cromwell programmer blunder: isSuccess was called on an incomplete RunStatus ($runStatus).")
    }
  }

  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] = {
    runStatus match {
      case successStatus: Succeeded => successStatus
        .eventList
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }
  }

  override def getTerminalMetadata(runStatus: RunStatus): Map[String, Any] = {
    runStatus match {
      case _: TerminalRunStatus => Map()
      case unknown => throw new RuntimeException(s"Attempt to get terminal metadata from non terminal status: $unknown")
    }
  }

  override val gcpBatchActor: ActorRef = backendSingletonActor

}




