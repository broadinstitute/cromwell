package cromwell.backend.google.batch.actors

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.batch.api.BatchApiRequestManager.BatchApiRequest
import cromwell.backend.google.batch.api.request.RequestHandler
import cromwell.backend.google.batch.api.BatchApiRequestManager
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.Mailbox
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive

object GcpBatchBackendSingletonActor {
  def props(qps: Int Refined Positive, requestWorkers: Int Refined Positive, serviceRegistryActor: ActorRef)(implicit
    requestHandler: RequestHandler
  ): Props =
    Props(
      new GcpBatchBackendSingletonActor(
        qps = qps,
        requestWorkers = requestWorkers,
        serviceRegistryActor = serviceRegistryActor
      )
    ).withDispatcher(BackendDispatcher)
}

final class GcpBatchBackendSingletonActor(
  qps: Int Refined Positive,
  requestWorkers: Int Refined Positive,
  serviceRegistryActor: ActorRef
)(implicit requestHandler: RequestHandler)
    extends Actor
    with ActorLogging {

  private val jesApiQueryManager = context.actorOf(
    BatchApiRequestManager.props(qps, requestWorkers, serviceRegistryActor).withMailbox(Mailbox.PriorityMailbox),
    "BatchQueryManager"
  )

  override def receive = {
    // Cromwell sends this message
    case abort: BackendSingletonActorAbortWorkflow =>
      // It seems that BatchAbortRequest is processed before this message, hence, we don't need to do anything else.
      // If it ever becomes necessary, we'll need to create link submitted jobs to its workflow id, which require
      // us to be cautious because batch deletes jobs instead of canceling them, hence, we should not delete jobs
      // that are on a final state.
      log.info(s"Cromwell requested to abort workflow ${abort.workflowId}")

    // TODO: Alex - is this still relevant?
    // jesApiQueryManager.forward(abort)

    case apiQuery: BatchApiRequest =>
      log.debug("Forwarding API query to PAPI request manager actor")
      jesApiQueryManager.forward(apiQuery)

    case other =>
      log.error(s"Unexpected message from {} to ${this.getClass.getSimpleName}: {}", sender().path.name, other)
  }
}
