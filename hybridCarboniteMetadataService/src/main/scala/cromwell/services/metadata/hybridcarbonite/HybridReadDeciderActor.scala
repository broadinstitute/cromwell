package cromwell.services.metadata.hybridcarbonite

import akka.actor.FSM.Normal
import akka.actor.{ActorRef, FSM, LoggingFSM, Props}
import cromwell.services.metadata.MetadataService.{MetadataReadAction, MetadataServiceResponse, QueryForWorkflowsMatchingParameters, WorkflowMetadataReadAction, WorkflowQueryFailure, WorkflowQuerySuccess}
import cromwell.services.metadata.hybridcarbonite.HybridReadDeciderActor._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.FailedMetadataResponse
import cromwell.services.metadata.{MetadataArchiveStatus, WorkflowQueryKey}

import scala.concurrent.{ExecutionContext, Future}

class HybridReadDeciderActor(classicMetadataServiceActor: ActorRef, carboniteMetadataServiceActor: ActorRef) extends LoggingFSM[HybridReadDeciderState, HybridReadDeciderData] {

  implicit val ec: ExecutionContext = context.dispatcher

  // TODO: [CARBONITE] Decide which actor to send the read request to
  // NOTE: the 'Future' return value was arbitrary to demonstrate the concept... feel free to refactor all this logic into something better
  //  like (eg) sending a message to the summary table followed by some appropriate forwarding action when the response comes back
  def decide(read: MetadataReadAction): Future[ActorRef] = Future.successful(classicMetadataServiceActor)


  when(Pending) {
    case Event(read: MetadataReadAction, NoData) =>
      val sndr = sender()

      read match {
        case _: WorkflowMetadataReadAction =>
          classicMetadataServiceActor ! QueryForWorkflowsMatchingParameters(Vector(WorkflowQueryKey.Id.name -> ""))
          goto(RequestingMetadataArchiveStatus) using WorkingData(sndr, read)

        case query: QueryForWorkflowsMatchingParameters =>
          classicMetadataServiceActor ! query
          goto(WaitingForMetadataResponse) using WorkingData(sndr, read)
      }
  }

  when(RequestingMetadataArchiveStatus) {
    case Event(WorkflowQuerySuccess(response, _), wd: WorkingData) =>
      if (response.results.size == 1 && response.results.headOption.exists(_.metadataArchiveStatus == MetadataArchiveStatus.Archived)) {
        carboniteMetadataServiceActor ! wd.request
      } else {
        classicMetadataServiceActor ! wd.request
      }
      goto(WaitingForMetadataResponse)

    case Event(WorkflowQueryFailure(reason), wd: WorkingData) =>
      log.error(reason, s"Failed to determine how to route ${wd.request.getClass.getSimpleName}. Falling back to classic.")
      classicMetadataServiceActor ! wd.request
      goto(WaitingForMetadataResponse)
  }

  when(WaitingForMetadataResponse) {
    case Event(response: MetadataServiceResponse , wd: WorkingData) =>
      wd.actor ! response
      stop(Normal)
  }

  whenUnhandled {
    case Event(msg, data) =>
      val errorMsg = s"Programmer Error: Unexpected message '$msg' sent to ${getClass.getSimpleName} in state $stateName with data $stateData from $sender()"
      log.error(errorMsg)
      data match {
        case NoData => stop(FSM.Failure(errorMsg))
        case WorkingData(actor, request) =>
          actor ! makeAppropriateFailureForRequest(errorMsg, request)

          stop(FSM.Failure(msg))
      }
  }

  def makeAppropriateFailureForRequest(msg: String, request: MetadataReadAction) = request match {
    case _: QueryForWorkflowsMatchingParameters => WorkflowQueryFailure(new Exception(msg))
    case _ => FailedMetadataResponse(request, new Exception(msg))
  }

}

object HybridReadDeciderActor {
  def props(classicMetadataServiceActor: ActorRef, carboniteMetadataServiceActor: ActorRef) =
    Props(new HybridReadDeciderActor(classicMetadataServiceActor, carboniteMetadataServiceActor))

  sealed trait HybridReadDeciderState
  case object Pending extends HybridReadDeciderState
  case object RequestingMetadataArchiveStatus extends HybridReadDeciderState
  case object WaitingForMetadataResponse extends HybridReadDeciderState

  sealed trait HybridReadDeciderData
  case object NoData extends HybridReadDeciderData
  final case class WorkingData(actor: ActorRef, request: MetadataReadAction) extends HybridReadDeciderData
}
