package cromwell.services.metadata.impl

import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core.WorkflowId
import cromwell.services.MetadataJsonResponse
import cromwell.services.metadata.MetadataService
import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, BuildWorkflowMetadataJsonAction, MetadataQueryResponse, MetadataServiceAction, MetadataServiceResponse, RootAndSubworkflowLabelsLookupResponse}
import cromwell.services.metadata.impl.ReadMetadataRegulatorActor.PropsMaker
import cromwell.services.metadata.impl.builder.MetadataBuilderActor

import scala.collection.mutable

class ReadMetadataRegulatorActor(metadataBuilderActorProps: PropsMaker, readMetadataWorkerProps: PropsMaker) extends Actor with ActorLogging {
  // This actor tracks all requests coming in from the API service and spins up new builders as needed to service them.
  // If the processing of an identical request is already in flight the requester will be added to a list of requesters
  // to notify when the response from the first request becomes available.
  // Note that if someone sends the same action twice, this actor must reply twice too.

  // Map from requests (MetadataServiceActions) to requesters.
  val apiRequests = new mutable.HashMap[MetadataServiceAction, List[ActorRef]]()
  // Map from ActorRefs of MetadataBuilderActors to requests. When a response comes back from a MetadataBuilderActor its
  // ActorRef is used as the lookup key in this Map. The result of that lookup yields the request which in turn is used
  // as the lookup key for requesters in the above Map.
  val builderRequests = new mutable.HashMap[ActorRef, MetadataServiceAction]()

  def createMetadataBuilderActor(workflowId: WorkflowId): ActorRef =
    context.actorOf(metadataBuilderActorProps().withDispatcher(ApiDispatcher), MetadataBuilderActor.uniqueActorName(workflowId.toString))

  def createReadMetadataWorker(): ActorRef =
    context.actorOf(readMetadataWorkerProps.apply().withDispatcher(ApiDispatcher), s"MetadataQueryWorker-${UUID.randomUUID()}")

  override def receive: Receive = {
    // This indirection via 'MetadataReadAction' lets the compiler make sure we cover all cases in the sealed trait:
    case action: BuildMetadataJsonAction =>
      action match {
        case singleWorkflowAction: BuildWorkflowMetadataJsonAction =>
          val currentRequesters = apiRequests.getOrElse(singleWorkflowAction, Nil)
          apiRequests.put(singleWorkflowAction, sender() +: currentRequesters)
          if (currentRequesters.isEmpty) {

            val builderActor = createMetadataBuilderActor(singleWorkflowAction.workflowId)
            builderRequests.put(builderActor, singleWorkflowAction)
            builderActor ! singleWorkflowAction
          }
        case crossWorkflowAction: MetadataService.QueryForWorkflowsMatchingParameters =>
          val currentRequesters = apiRequests.getOrElse(crossWorkflowAction, Nil)
          apiRequests.put(crossWorkflowAction, sender() +: currentRequesters)
          if (currentRequesters.isEmpty) {
            val readMetadataActor = createReadMetadataWorker()
            builderRequests.put(readMetadataActor, crossWorkflowAction)
            readMetadataActor ! crossWorkflowAction
          }
      }
    case serviceResponse: MetadataServiceResponse =>
      serviceResponse match {
        case response @ (_: MetadataJsonResponse | _: MetadataQueryResponse | _: RootAndSubworkflowLabelsLookupResponse) =>
          handleResponseFromMetadataWorker(response)
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
  type PropsMaker = () => Props

  def props(singleWorkflowMetadataBuilderProps: PropsMaker, summarySearcherProps: PropsMaker): Props = {
    Props(new ReadMetadataRegulatorActor(singleWorkflowMetadataBuilderProps, summarySearcherProps))
  }
}
