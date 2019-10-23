package cromwell.services.metadata.hybridcarbonite

import akka.actor.{ActorRef, FSM, LoggingFSM, Props}
import cromwell.services.FailedMetadataResponse
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.hybridcarbonite.HybridReadDeciderActor._
import cromwell.services.metadata.{MetadataArchiveStatus, WorkflowQueryKey}

import scala.concurrent.ExecutionContext

class HybridReadDeciderActor(classicMetadataServiceActor: ActorRef, carboniteMetadataServiceActor: ActorRef) extends LoggingFSM[HybridReadDeciderState, HybridReadDeciderData] {

  startWith(Pending, NoData)

  implicit val ec: ExecutionContext = context.dispatcher

  when(Pending) {
    case Event(action: MetadataReadAction, NoData) if action.requiresOnlyClassicMetadata =>
      classicMetadataServiceActor ! action
      goto(WaitingForMetadataResponse) using WorkingData(sender(), action)
    case Event(action: WorkflowMetadataReadAction, NoData) =>
      classicMetadataServiceActor ! QueryForWorkflowsMatchingParameters(Vector(WorkflowQueryKey.Id.name -> action.workflowId.toString))
      goto(RequestingMetadataArchiveStatus) using WorkingData(sender(), action)
  }

  when(RequestingMetadataArchiveStatus) {
    case Event(s: WorkflowQuerySuccess, wd: WorkingData) if s.hasMultipleSummaryRows =>
      val errorMsg = s"Programmer Error: Found multiple summary rows looking up metadata archive status for ${wd.request}: ${s.response}"
      wd.requester ! makeAppropriateFailureForRequest(errorMsg, wd.request)
      stop(FSM.Failure(errorMsg))
    case Event(s: WorkflowQuerySuccess, wd: WorkingData) if s.isCarbonited =>
      carboniteMetadataServiceActor ! wd.request
      goto(WaitingForMetadataResponse)
    case Event(_: WorkflowQuerySuccess, wd: WorkingData) => // this is an uncarbonited workflow
      classicMetadataServiceActor ! wd.request
      goto(WaitingForMetadataResponse)
    case Event(WorkflowQueryFailure(reason), wd: WorkingData) =>
      log.error(reason, s"Programmer Error: Failed to determine how to route ${wd.request.getClass.getSimpleName}. Falling back to classic.")
      classicMetadataServiceActor ! wd.request
      goto(WaitingForMetadataResponse)
  }

  when(WaitingForMetadataResponse) {
    case Event(response: MetadataServiceResponse, wd: WorkingData) =>
      wd.requester ! response
      stop(FSM.Normal)
  }

  whenUnhandled {
    case Event(akkaMessage, data) =>
      val logMessage = s"Programmer Error: Unexpected message '$akkaMessage' sent to ${getClass.getSimpleName} in state $stateName with data $stateData from $sender()"
      log.error(logMessage)
      data match {
        case NoData =>
          stop(FSM.Failure(logMessage))
        case WorkingData(actor, request) =>
          actor ! makeAppropriateFailureForRequest(logMessage, request)
          stop(FSM.Failure(logMessage))
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
  final case class WorkingData(requester: ActorRef, request: MetadataReadAction) extends HybridReadDeciderData

  implicit class EnhancedWorkflowQuerySuccess(val success: WorkflowQuerySuccess) extends AnyVal {
    def hasMultipleSummaryRows: Boolean = success.response.results.size > 1
    def isCarbonited: Boolean = success.response.results.headOption.exists(_.metadataArchiveStatus == MetadataArchiveStatus.Archived)
  }

  implicit class EnhancedMetadataReadAction(val action: MetadataReadAction) extends AnyVal {
    def requiresOnlyClassicMetadata: Boolean = action match {
      case _: GetLabels | _: GetRootAndSubworkflowLabels | _: GetStatus | _: QueryForWorkflowsMatchingParameters => true
      case _ => false // GetLogs, WorkflowOutputs, GetMetadataAction
    }
  }
}
