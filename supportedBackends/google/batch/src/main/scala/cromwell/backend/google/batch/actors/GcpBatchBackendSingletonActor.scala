package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.batch.api.BatchApiRequestManager
import cromwell.backend.google.batch.api.BatchApiRequestManager.BatchApiRequest
import cromwell.backend.google.batch.api.request.{BatchRequestExecutor, RequestHandler}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.Mailbox
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive

object GcpBatchBackendSingletonActor {
  def props(qps: Int Refined Positive,
            requestWorkers: Int Refined Positive,
            serviceRegistryActor: ActorRef,
            batchRequestExecutor: BatchRequestExecutor
  )(implicit
    requestHandler: RequestHandler
  ): Props =
    Props(
      new GcpBatchBackendSingletonActor(
        qps = qps,
        requestWorkers = requestWorkers,
        serviceRegistryActor = serviceRegistryActor,
        batchRequestExecutor = batchRequestExecutor
      )
    ).withDispatcher(BackendDispatcher)
}

final class GcpBatchBackendSingletonActor(
  qps: Int Refined Positive,
  requestWorkers: Int Refined Positive,
  serviceRegistryActor: ActorRef,
  batchRequestExecutor: BatchRequestExecutor
)(implicit requestHandler: RequestHandler)
    extends Actor
    with ActorLogging {

  private val batchApiRequestManager = context.actorOf(
    BatchApiRequestManager
      .props(qps = qps, requestWorkers = requestWorkers, serviceRegistryActor, batchRequestExecutor)
      .withMailbox(Mailbox.PriorityMailbox),
    "BatchApiRequestManager"
  )

  override def receive = {
    case _: BackendSingletonActorAbortWorkflow =>
      // Cromwell sends this message to abort a job, but "GcpBatchAsyncBackendJobExecutionActor" implements
      // the "tryAbort" method which makes this message irrelevant.
      //
      // If it ever becomes necessary, we'll need to link submitted jobs to its workflow id, which require
      // us to be cautious because batch deletes jobs instead of canceling them, hence, we should not delete jobs
      // that are on a final state.
      //
      // batchApiRequestManager.forward(abort)
      ()

    case apiQuery: BatchApiRequest =>
      log.debug("Forwarding API query to Batch request manager actor")
      batchApiRequestManager.forward(apiQuery)

    case other =>
      log.error(s"Unexpected message from {} to ${this.getClass.getSimpleName}: {}", sender().path.name, other)
  }
}
