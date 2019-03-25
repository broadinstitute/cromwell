package cromwell.services.metadata.impl


import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.impl.MetadataSummaryRefreshActor._

import scala.util.{Failure, Success}

/**
  * This looks for workflows whose metadata summaries are in need of refreshing and refreshes those summaries.
  * Despite its package location this is not actually a service, but a type of actor spawned at the behest of
  * the MetadataServiceActor.
  *
  */

object MetadataSummaryRefreshActor {
  sealed trait MetadataSummaryActorMessage
  final case class SummarizeMetadata(limit: Int, respondTo: ActorRef) extends MetadataSummaryActorMessage
  case object MetadataSummarySuccess extends MetadataSummaryActorMessage
  final case class MetadataSummaryFailure(t: Throwable) extends MetadataSummaryActorMessage

  def props() = Props(new MetadataSummaryRefreshActor()).withDispatcher(ServiceDispatcher)

  sealed trait SummaryRefreshState
  case object WaitingForRequest extends SummaryRefreshState
  case object SummarizingMetadata extends SummaryRefreshState
  case object MetadataSummaryComplete extends SummaryRefreshState

  case object SummaryRefreshData
}

class MetadataSummaryRefreshActor()
  extends LoggingFSM[SummaryRefreshState, SummaryRefreshData.type]
    with MetadataDatabaseAccess with MetadataServicesStore {

  implicit val ec = context.dispatcher

  startWith(WaitingForRequest, SummaryRefreshData)

  when (WaitingForRequest) {
    case Event(SummarizeMetadata(limit, respondTo), _) =>
      refreshWorkflowMetadataSummaries(limit) onComplete {
        case Success(summaryResult) =>
          log.debug(
            "Summarized metadata: increasing = {}, decreasing = {}",
            summaryResult.increasingResult,
            summaryResult.decreasingResult
          )
          respondTo ! MetadataSummarySuccess
          self ! MetadataSummaryComplete
        case Failure(t) =>
          log.error(t, "Failed to summarize metadata")
          respondTo ! MetadataSummaryFailure(t)
          self ! MetadataSummaryComplete
      }
      goto(SummarizingMetadata)
  }

  when (SummarizingMetadata) {
    case Event(MetadataSummaryComplete, _) =>
      goto(WaitingForRequest) using SummaryRefreshData
  }

  whenUnhandled {
    case Event(wut, _) =>
      log.warning("Unrecognized or unexpected message while in state '{}': {}", stateName, wut)
      stay()
  }
}
