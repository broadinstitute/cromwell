package cromwell.services.metadata.impl


import akka.actor.{ActorRef, LoggingFSM, Props}
import com.typesafe.config.ConfigFactory
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
  final case class SummarizeMetadata(respondTo: ActorRef) extends MetadataSummaryActorMessage
  case object MetadataSummarySuccess extends MetadataSummaryActorMessage
  final case class MetadataSummaryFailure(t: Throwable) extends MetadataSummaryActorMessage

  final case class MetadataSummaryCompleteMessage(newData: SummaryRefreshData) extends MetadataSummaryActorMessage

  def props() = Props(new MetadataSummaryRefreshActor()).withDispatcher(ServiceDispatcher)

  sealed trait SummaryRefreshState
  case object WaitingForRequest extends SummaryRefreshState
  case object SummarizingMetadata extends SummaryRefreshState

  case class SummaryRefreshData(mostRecentlyLoggedTotal: Option[Long])
}

class MetadataSummaryRefreshActor()
  extends LoggingFSM[SummaryRefreshState, SummaryRefreshData]
    with MetadataDatabaseAccess with MetadataServicesStore {

  val config = ConfigFactory.load
  implicit val ec = context.dispatcher

  startWith(WaitingForRequest, SummaryRefreshData(None))

  when (WaitingForRequest) {
    case Event(SummarizeMetadata(respondTo), oldData @ SummaryRefreshData(mostRecentlyLoggedTotal)) =>
      refreshWorkflowMetadataSummaries() onComplete {
        case Success(newMaximumSummaryId) =>
          val newData = if(mostRecentlyLoggedTotal.contains(newMaximumSummaryId)) {
            oldData
          } else {
            log.info(s"Metadata summarizer has now reached: $newMaximumSummaryId")
            SummaryRefreshData(Some(newMaximumSummaryId))
          }

          respondTo ! MetadataSummarySuccess
          self ! MetadataSummaryCompleteMessage(newData)
        case Failure(t) =>
          log.error(t, "Failed to summarize metadata")
          respondTo ! MetadataSummaryFailure(t)
          self ! MetadataSummaryCompleteMessage(oldData)
      }
      goto(SummarizingMetadata)
  }

  when (SummarizingMetadata) {
    case Event(MetadataSummaryCompleteMessage(newData), _) =>
      goto(WaitingForRequest) using newData
  }

  whenUnhandled {
    case Event(wut, _) =>
      log.warning("Unrecognized or unexpected message while in state '{}': {}", stateName, wut)
      stay()
  }
}
