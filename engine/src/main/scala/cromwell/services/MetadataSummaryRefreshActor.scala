package cromwell.services

import java.time.OffsetDateTime

import akka.actor.{LoggingFSM, Props}
import com.typesafe.config.ConfigFactory
import cromwell.engine.db.DataAccess
import cromwell.services.MetadataSummaryRefreshActor._

import scala.util.{Failure, Success}


/**
  * This looks for workflows whose metadata summaries are in need of refreshing and refreshes those summaries.
  * Despite its package location this is not actually a service, but a type of actor spawned at the behest of
  * the EngineMetadataServiceActor.
  *
  */

object MetadataSummaryRefreshActor {
  sealed trait MetadataSummaryActorMessage
  case object SummarizeMetadata
  case object MetadataSummarySuccess
  case class MetadataSummaryFailure(t: Throwable)

  def props(startMetadataTimestamp: Option[OffsetDateTime]) = Props(new MetadataSummaryRefreshActor(startMetadataTimestamp))

  trait SummaryRefreshState
  case object WaitingForRequest extends SummaryRefreshState
  case object SummarizingMetadata extends SummaryRefreshState
  private case class MetadataSummaryComplete(startMetadataId: Long) extends SummaryRefreshState

  case class SummaryRefreshData(startMetadataId: Long)
}


class MetadataSummaryRefreshActor(startMetadataTimestamp: Option[OffsetDateTime]) extends LoggingFSM[SummaryRefreshState, SummaryRefreshData] {

  val dataAccess = DataAccess.globalDataAccess
  val config = ConfigFactory.load
  implicit val ec = context.dispatcher

  startWith(WaitingForRequest, SummaryRefreshData(startMetadataId = 0L))

  when (WaitingForRequest) {
    case (Event(SummarizeMetadata, data)) =>
      val sndr = sender()
      val startMetadataId = data.startMetadataId
      dataAccess.refreshWorkflowMetadataSummaries(startMetadataId, startMetadataTimestamp) onComplete {
        case Success(id) =>
          sndr ! MetadataSummarySuccess
          self ! MetadataSummaryComplete(startMetadataId = id + 1)
        case Failure(t) =>
          log.error(t, "Failed to summarize metadata starting from index {}", startMetadataId)
          sndr ! MetadataSummaryFailure(t)
          self ! MetadataSummaryComplete(startMetadataId = startMetadataId)
      }
      goto(SummarizingMetadata)
  }

  when (SummarizingMetadata) {
    case Event(MetadataSummaryComplete(startMetadataId), _) =>
      goto(WaitingForRequest) using SummaryRefreshData(startMetadataId)
  }

  whenUnhandled {
    case Event(wut, _) =>
      log.warning("Unrecognized or unexpected message while in state '{}': {}", stateName, wut)
      stay()
  }
}
