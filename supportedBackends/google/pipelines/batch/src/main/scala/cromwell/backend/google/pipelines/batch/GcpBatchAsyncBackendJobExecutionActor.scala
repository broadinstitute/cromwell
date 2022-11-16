package cromwell.backend.google.pipelines.batch

//import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.backend.google.pipelines.common

import scala.concurrent.duration._
import cromwell.backend._
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.google.pipelines.common.Run
//import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.async.ExecutionHandle
//import cromwell.backend.google.pipelines.common.Run
//import cromwell.core._

import scala.concurrent.Future
//import scala.concurrent.{Future, Promise}
//import cromwell.backend.async.ExecutionHandle

import akka.actor.ActorRef
//import akka.actor.{Actor, ActorRef}
//import cromwell.backend.google.pipelines.batch.GcpBatchRunCreationClient

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams) extends BackendJobLifecycleActor with StandardAsyncExecutionActor with GcpBatchRunCreationClient {
  /** The type of the run info when a job is started. */
  //override type StandardAsyncRunInfo = this.type
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = this.type
  //override type StandardAsyncRunState = RunStatus


  //import GcpBatchAsyncBackendJobExecutionActor._

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  override def statusEquivalentTo(thiz: GcpBatchAsyncBackendJobExecutionActor
    .this.type)(that: GcpBatchAsyncBackendJobExecutionActor
    .this.type): Boolean = true

  /**
    * Returns true when a job is complete, either successfully or unsuccessfully.
    *
    * @param runStatus The run status.
    * @return True if the job has completed.
    */
  override def isTerminal(runStatus: GcpBatchAsyncBackendJobExecutionActor
    .this.type): Boolean = true

  override def dockerImageUsed: Option[String] = Option("test")

  //lazy val batchJob: GcpBatchJob = GcpBatchJob(jobDescriptor)
  val backendSingletonActor: ActorRef

  // Primary entry point for cromwell to run GCP Batch job
  override def executeAsync(): Future[ExecutionHandle] = {


    //val runId = backendSingletonActor ! runPipeline()
    val runId = runPipeline()
    //backendSingletonActor.notify()
    //PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None)
    //ExecutionHandle
    //val runId: StandardAsyncJob
    PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None)

  }

  //yield PendingExecutionHandle(jobDescriptor = jobDescriptor, previousState = None, job = StandardAsyncJob("id"))

  //override def pollBackOff: SimpleExponentialBackoff = ???

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(1
    .second, 5
    .minutes, 1.1)

  //override def executeOrRecoverBackOff: SimpleExponentialBackoff = ???

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 3
      .second, maxInterval = 20
      .second, multiplier = 1.1)


}

//case class PendingExecutionHandle(jobDescriptor: BackendJobDescriptor, previousState: None.type, job: StandardAsyncJob)
