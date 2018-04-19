package cromwell.backend.google.pipelines.common

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.pipelines.common.api.{PipelinesApiBatchHandler, PipelinesApiRequestManager}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIApiRequest
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.Mailbox
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive

final case class PipelinesApiBackendSingletonActor(qps: Int Refined Positive, requestWorkers: Int Refined Positive, serviceRegistryActor: ActorRef)
                                                  (implicit batchHandler: PipelinesApiBatchHandler) extends Actor with ActorLogging {

  val jesApiQueryManager = context.actorOf(PipelinesApiRequestManager.props(qps, requestWorkers, serviceRegistryActor).withMailbox(Mailbox.PriorityMailbox), "PAPIQueryManager")

  override def receive = {
    case abort: BackendSingletonActorAbortWorkflow => jesApiQueryManager.forward(abort)
    case apiQuery: PAPIApiRequest =>
      log.debug("Forwarding API query to JES API query manager actor")
      jesApiQueryManager.forward(apiQuery)
  }
}

object PipelinesApiBackendSingletonActor {
  def props[O](qps: Int Refined Positive, requestWorkers: Int Refined Positive, serviceRegistryActor: ActorRef)
              (implicit batchHandler: PipelinesApiBatchHandler): Props = Props(PipelinesApiBackendSingletonActor(qps, requestWorkers, serviceRegistryActor)).withDispatcher(BackendDispatcher)
}
