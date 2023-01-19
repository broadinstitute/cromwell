package cromwell.backend.google.pipelines.batch

import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.backend._
import cromwell.backend.google.pipelines.batch.RunStatus.{DeletionInProgress, Failed, StateUnspecified, Unrecognized}
//import cromwell.backend.google.pipelines.common.WorkflowOptionKeys
import cromwell.core.{ExecutionEvent, WorkflowId}
//import wom.callable.RuntimeEnvironment
//import cromwell.core.{ExecutionEvent, WorkflowId}
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.async.ExecutionHandle
import akka.actor.ActorRef
import akka.pattern.AskSupport
import cromwell.services.instrumentation.CromwellInstrumentation

import java.util.UUID
import scala.concurrent.Future
import scala.concurrent.duration._
import GcpBatchBackendSingletonActor._
import cromwell.backend.google.pipelines.batch.RunStatus.{Running,Succeeded, TerminalRunStatus}

import com.google.cloud.batch.v1.JobStatus

//import scala.util.Success
//import scala.util.{Failure, Success, Try}
//import cromwell.core.ExecutionEvent

object GcpBatchAsyncBackendJobExecutionActor {

  type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, RunStatus]

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams) extends BackendJobLifecycleActor with StandardAsyncExecutionActor with AskSupport with GcpBatchJobCachingActorHelper with GcpBatchStatusRequestClient with CromwellInstrumentation {

  import GcpBatchAsyncBackendJobExecutionActor._

  lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id
  //val requestFactory: PipelinesApiRequestFactory = initializationData.genomicsRequestFactory

  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = RunStatus

  // temporary until GCP Batch client can generate random job IDs
  private val jobTemp = "job-" + java.util.UUID.randomUUID.toString

  //override def receive: Receive = pollingActorClientReceive orElse runCreationClientReceive orElse super.receive

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  //override def dockerImageUsed: Option[String] = Option("test")

  private lazy val jobDockerImage = jobDescriptor.maybeCallCachingEligible.dockerHash.getOrElse(runtimeAttributes.dockerImage)
  override def dockerImageUsed: Option[String] = Option(jobDockerImage)

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

    val batchTest = BatchRequest(projectId="batch-testing-350715", region="us-central1", jobName=jobTemp, dockerImage=jobDockerImage)

    val runBatchResponse = for {
      _ <- uploadScriptFile()
      _ = backendSingletonActor ! batchTest
      runId = StandardAsyncJob(UUID.randomUUID().toString)  //temp to test
    }
    yield runId

    runBatchResponse map {runId => PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None)}

  }

  override def reconnectAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = {
    log.info("reconnect async runs") // in for debugging remove later
    val handle = PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](jobDescriptor, jobId, Option(Run(jobId)), previousState = None)
    Future.successful(handle)
  }

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(1
    .second, 5
    .minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 3
      .second, maxInterval = 20
      .second, multiplier = 1.1)

  override def pollStatusAsync(handle: GcpBatchPendingExecutionHandle): Future[RunStatus] = {

     //println(batchAttributes.project)
     val gcpBatchPoll = new GcpBatchJobGetRequest
     val result = gcpBatchPoll.GetJob(jobTemp)
     //val temp = result.toString //matches for string
     //val batchRunStatus = RunStatus.fromJobStatus(status=result)
     //val eventList: Seq[ExecutionEvent] = Seq(ExecutionEvent.toString)
     val jobStatus = result.getStatus.getState

    //https://github.com/broadinstitute/cromwell/blob/328a0fe0aa307ee981b00e4af6b397b61a9fbe9e/engine/src/main/scala/cromwell/engine/workflow/lifecycle/execution/SubWorkflowExecutionActor.scala
     jobStatus match {
       case JobStatus.State.QUEUED =>
         log.info("job queued")
         Future.successful(Running)
       case JobStatus.State.SCHEDULED =>
         log.info("job scheduled")
         Future.successful(Running)
       case JobStatus.State.RUNNING =>
         log.info("job running")
         Future.successful(Running)
       case JobStatus.State.SUCCEEDED =>
         log.info("job scheduled")
         Future.successful(Succeeded(List(ExecutionEvent("complete in GCP Batch")))) //update to more specific
       case JobStatus.State.FAILED =>
         log.info("job failed")
         Future.successful(Failed)
       case JobStatus.State.DELETION_IN_PROGRESS =>
         log.info("deletion in progress")
         Future.successful(DeletionInProgress)
       case JobStatus.State.STATE_UNSPECIFIED =>
         log.info("state unspecified")
         Future.successful(StateUnspecified)
       case JobStatus.State.UNRECOGNIZED =>
         log.info("state unrecognized")
         Future.successful(Unrecognized)
       //case _ =>
       //  log.info("job status not matched")
       //  Future.successful(Running)

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
        println(f"isTerminal match TerminalRunStatus running with status $tempTermStatus")
        true
      case other =>
        println(f"isTerminal match _ running with status $other")
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



  /*
  override lazy val runtimeEnvironment: RuntimeEnvironment = {
    RuntimeEnvironmentBuilder(jobDescriptor
      .runtimeAttributes, PipelinesApiWorkingDisk
      .MountPoint, PipelinesApiWorkingDisk
      .MountPoint)(standardParams
      .minimumRuntimeSettings)
  }
  */

}




