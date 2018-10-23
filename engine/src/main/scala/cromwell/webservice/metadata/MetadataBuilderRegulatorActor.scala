package cromwell.webservice.metadata

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.effect.IO
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.services.metadata.MetadataService.{MetadataServiceAction, StreamedGetMetadataQueryAction, StreamedGetSingleWorkflowMetadataAction}
import cromwell.webservice.metadata.MetadataBuilderActor.{BuiltMetadataResponse, FailedMetadataResponse, MetadataBuilderActorResponse}
import spray.json.JsObject

import scala.collection.mutable
import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success}

class MetadataBuilderRegulatorActor(serviceRegistryActor: ActorRef, streamTimeout: FiniteDuration) extends Actor with ActorLogging {
  // This actor tracks all requests coming in from the API service and spins up new builders as needed to service them.
  // If the processing of an identical request is already in flight the requester will be added to a set of requesters
  // to notify when the response from the first request becomes available.

  // Map from requests (MetadataServiceActions) to requesters.
  val apiRequests = new mutable.HashMap[MetadataServiceAction, Set[ActorRef]]()
  // Map from ActorRefs of MetadataBuilderActors to requests. When a response comes back from a MetadataBuilderActor its
  // ActorRef is used as the lookup key in this Map. The result of that lookup yields the request which in turn is used
  // as the lookup key for requesters in the above Map.
  val builderRequests = new mutable.HashMap[ActorRef, MetadataServiceAction]()

  private val streamMetadataBuilder = new StreamMetadataBuilder(streamTimeout)

  implicit lazy val ec = context.system.dispatchers.lookup(ApiDispatcher)

  override def receive: Receive = {
    case action@StreamedGetSingleWorkflowMetadataAction(workflowId, includeKeysOption, excludeKeysOption, expandSubWorkflows) =>
      handleStreamAction(action, streamMetadataBuilder.workflowMetadataQuery(workflowId, includeKeysOption, excludeKeysOption, expandSubWorkflows))
    case action: StreamedGetMetadataQueryAction =>
      handleStreamAction(action, streamMetadataBuilder.workflowMetadataQuery(action.key))
    case action: MetadataServiceAction =>
      handleAction(action) {
        val metadataBuilderActor = context.actorOf(
          MetadataBuilderActor.props(serviceRegistryActor).withDispatcher(ApiDispatcher), MetadataBuilderActor.uniqueActorName)
        builderRequests.put(metadataBuilderActor, action)
        metadataBuilderActor ! action
      }
    case (action: MetadataServiceAction, response: MetadataBuilderActorResponse) =>
      processResponse(action, response)
    case response: MetadataBuilderActorResponse =>
      val sndr = sender()
      builderRequests.get(sndr) match {
        case Some(action) => processResponse(action, response)
        case None =>
          // unpossible: this actor should know about all the child MetadataBuilderActors it has begotten
          log.error(s"MetadataBuilderRegulatorActor unpossible error: unrecognized sender $sndr")
      }
  }

  private def processResponse(action: MetadataServiceAction, response: MetadataBuilderActorResponse): Unit = {
    apiRequests.get(action) match {
      case Some(requesters) =>
        apiRequests.remove(action)
        requesters foreach {
          _ ! response
        }
      case None =>
        // unpossible: there had to have been a request that corresponded to this response
        log.error(s"MetadataBuilderRegulatorActor unpossible error: no requesters found for action: $action")
    }
    builderRequests.remove(sender())
    ()
  }

  private def handleAction(action: MetadataServiceAction)(onEmpty: => Unit): Unit = {
    val currentRequesters = apiRequests.getOrElse(action, Set.empty)
    apiRequests.put(action, currentRequesters + sender())
    if (currentRequesters.isEmpty) onEmpty
  }

  private def handleStreamAction(action: MetadataServiceAction, ioJson: IO[JsObject]) = {
    handleAction(action) {
      ioJson.unsafeToFuture().map(BuiltMetadataResponse.apply) onComplete {
        case Success(v) =>
          self.tell((action, v), ActorRef.noSender)
        case Failure(f) =>
          self.tell((action, FailedMetadataResponse(f)), ActorRef.noSender)
      }
    }
  }
}

object MetadataBuilderRegulatorActor {
  def props(serviceRegistryActor: ActorRef, streamTimeout: FiniteDuration): Props = {
    Props(new MetadataBuilderRegulatorActor(serviceRegistryActor, streamTimeout))
  }
}
