package cromwell.services.metadata.impl

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.services.metadata.MetadataService.MetadataServiceAction
import cromwell.services.metadata.impl.ReadMetadataRegulatorActor.ReadMetadataWorkerMaker
import cromwell.services.metadata.impl.builder.MetadataBuilderActor
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.MetadataBuilderActorResponse

import scala.collection.mutable

class ReadMetadataRegulatorActor(readMetadataWorkerMaker: ReadMetadataWorkerMaker) extends Actor with ActorLogging {
  // This actor tracks all requests coming in from the API service and spins up new builders as needed to service them.
  // If the processing of an identical request is already in flight the requester will be added to a set of requesters
  // to notify when the response from the first request becomes available.

  // Map from requests (MetadataServiceActions) to requesters.
  val apiRequests = new mutable.HashMap[MetadataServiceAction, Set[ActorRef]]()
  // Map from ActorRefs of MetadataBuilderActors to requests. When a response comes back from a MetadataBuilderActor its
  // ActorRef is used as the lookup key in this Map. The result of that lookup yields the request which in turn is used
  // as the lookup key for requesters in the above Map.
  val builderRequests = new mutable.HashMap[ActorRef, MetadataServiceAction]()

  override def receive: Receive = {
    case action: MetadataServiceAction =>
      val currentRequesters = apiRequests.getOrElse(action, Set.empty)
      apiRequests.put(action, currentRequesters + sender())
      if (currentRequesters.isEmpty) {
        val readMetadataActor = context.actorOf(readMetadataWorkerMaker.apply().withDispatcher(ApiDispatcher), MetadataBuilderActor.uniqueActorName)
        builderRequests.put(readMetadataActor, action)
        readMetadataActor ! action
      }
    case response: MetadataBuilderActorResponse =>
      val sndr = sender()
      builderRequests.get(sndr) match {
        case Some(action) =>
          apiRequests.get(action) match {
            case Some(requesters) =>
              apiRequests.remove(action)
              requesters foreach { _ ! response}
            case None =>
              // unpossible: there had to have been a request that corresponded to this response
              log.error(s"MetadataBuilderRegulatorActor unpossible error: no requesters found for action: $action")
          }
          builderRequests.remove(sndr)
          ()
        case None =>
          // unpossible: this actor should know about all the child MetadataBuilderActors it has begotten
          log.error(s"MetadataBuilderRegulatorActor unpossible error: unrecognized sender $sndr")
      }
  }
}

object ReadMetadataRegulatorActor {

  type ReadMetadataWorkerMaker = () => Props

  def props(readMetadataWorkerMaker: ReadMetadataWorkerMaker): Props = {
    Props(new ReadMetadataRegulatorActor(readMetadataWorkerMaker))
  }
}
