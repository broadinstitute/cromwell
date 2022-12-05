package cromwell.backend.google.pipelines.batch

import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration._
import cromwell.backend._
import cromwell.backend.google.pipelines.common.Run
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.async.ExecutionHandle
import akka.actor.ActorRef
import akka.stream.ActorMaterializer

import java.util.UUID
import scala.concurrent.Future

object GcpBatchAsyncBackendJobExecutionActor {

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams) extends BackendJobLifecycleActor with StandardAsyncExecutionActor with GcpBatchRunCreationClient {

  //testing actor creation.  may not be needed
  implicit val actorSystem = context.system
  implicit val materializer = ActorMaterializer()

  /** The type of the run info when a job is started. */
  //override type StandardAsyncRunInfo = this.type
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = this.type
  //override type StandardAsyncRunState = RunStatus


  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  override def statusEquivalentTo(thiz: GcpBatchAsyncBackendJobExecutionActor
    .this.type)(that: GcpBatchAsyncBackendJobExecutionActor
    .this.type): Boolean = true

  override def isTerminal(runStatus: GcpBatchAsyncBackendJobExecutionActor
    .this.type): Boolean = true

  override def dockerImageUsed: Option[String] = Option("test")

  val backendSingletonActor: ActorRef = standardParams.backendSingletonActorOption.getOrElse(throw new RuntimeException("GCP Batch actor cannot exist without its backend singleton 2"))

  //lazy val batchJob: GcpBatchJob = GcpBatchJob(jobDescriptor)
  override val gcpBatchApiActor: ActorRef = standardParams.backendSingletonActorOption.getOrElse(
    throw new RuntimeException("Batch Backend actor cannot exist without the Batch backend singleton actor"))

  // Primary entry point for cromwell to run GCP Batch job
  //override def executeAsync(): Future[ExecutionHandle] = {
  override def executeAsync(): Future[ExecutionHandle] = {

      //backendSingletonActor ! runPipeline()
      runPipeline()
      val runId = StandardAsyncJob(UUID.randomUUID().toString)  //temp to test

      Future.successful(PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None))

      //ExecutionHandle[]
  }

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(1
    .second, 5
    .minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 3
      .second, maxInterval = 20
      .second, multiplier = 1.1)


}
