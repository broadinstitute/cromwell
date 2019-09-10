package cromwell.services.metadata.impl

import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.services.metadata.MetadataService
import cromwell.services.metadata.MetadataService.{MetadataQueryResponse, MetadataReadAction, MetadataServiceAction, MetadataServiceResponse, WorkflowMetadataReadAction}
import cromwell.services.metadata.impl.ReadMetadataRegulatorActor.ReadMetadataWorkerMaker
import cromwell.services.metadata.impl.builder.MetadataBuilderActor
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.MetadataBuilderActorResponse

import scala.collection.mutable

class ReadMetadataRegulatorActor(metadataBuilderActorProps: ReadMetadataWorkerMaker, readMetadataWorkerProps: ReadMetadataWorkerMaker) extends Actor with ActorLogging {
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
    // This indirection via 'MetadataReadAction' lets the compiler make sure we cover all cases in the sealed trait:
    case action: MetadataReadAction =>
      action match {
        case singleWorkflowAction: WorkflowMetadataReadAction =>
          val currentRequesters = apiRequests.getOrElse(singleWorkflowAction, Set.empty)
          apiRequests.put(singleWorkflowAction, currentRequesters + sender())
          if (currentRequesters.isEmpty) {

            val builderActor = context.actorOf(metadataBuilderActorProps().withDispatcher(ApiDispatcher), MetadataBuilderActor.uniqueActorName(singleWorkflowAction.workflowId.toString))
            builderRequests.put(builderActor, singleWorkflowAction)
            builderActor ! singleWorkflowAction
          }
        case crossWorkflowAction: MetadataService.QueryForWorkflowsMatchingParameters =>
          val currentRequesters = apiRequests.getOrElse(crossWorkflowAction, Set.empty)
          apiRequests.put(crossWorkflowAction, currentRequesters + sender())
          if (currentRequesters.isEmpty) {
            val readMetadataActor = context.actorOf(readMetadataWorkerProps.apply().withDispatcher(ApiDispatcher), s"MetadataQueryWorker-${UUID.randomUUID()}")
            builderRequests.put(readMetadataActor, crossWorkflowAction)
            readMetadataActor ! crossWorkflowAction
          }
      }
    case serviceResponse: MetadataServiceResponse =>
      serviceResponse match {
        case response: MetadataBuilderActorResponse => handleResponseFromMetadataWorker(response)
        case response: MetadataQueryResponse => handleResponseFromMetadataWorker(response)
      }
    case other => log.error(s"Programmer Error: Unexpected message $other received from $sender")
  }

  def handleResponseFromMetadataWorker(response: Any): Unit = {
    val sndr = sender()
    builderRequests.get(sndr) match {
      case Some(action) =>
        apiRequests.get(action) match {
          case Some(requesters) =>
            apiRequests.remove(action)
            requesters foreach { _ ! response}
          case None =>
            // unpossible: there had to have been a request that corresponded to this response
            log.error(s"Programmer Error: MetadataBuilderRegulatorActor has no registered requesters found for action: $action")
        }
        builderRequests.remove(sndr)
        ()
      case None =>
        // unpossible: this actor should know about all the child MetadataBuilderActors it has begotten
        log.error(s"Programmer Error: MetadataBuilderRegulatorActor received a metadata response from an unrecognized sender $sndr")
    }
  }
}

object ReadMetadataRegulatorActor {

  type ReadMetadataWorkerMaker = () => Props

  def props(metadataBuilderActorProps: ReadMetadataWorkerMaker, readMetadataWorkerProps: ReadMetadataWorkerMaker): Props = {
    Props(new ReadMetadataRegulatorActor(metadataBuilderActorProps, readMetadataWorkerProps))
  }
}
