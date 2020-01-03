package cromwell.services.metadata.impl


import akka.actor.{ActorRef, LoggingFSM, Props}
import cats.data.NonEmptyList
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
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

  def props(serviceRegistryActor: ActorRef) = Props(new MetadataSummaryRefreshActor(serviceRegistryActor)).withDispatcher(ServiceDispatcher)

  sealed trait SummaryRefreshState
  case object WaitingForRequest extends SummaryRefreshState
  case object SummarizingMetadata extends SummaryRefreshState

  sealed trait SummaryRefreshData
  case object EmptySummaryRefreshData extends SummaryRefreshData
  final case class PreviousMaximumMetadataEntryId(value: Long) extends SummaryRefreshData

  // Internal message to self
  case class MetadataSummaryComplete(nextState: SummaryRefreshData)
}

class MetadataSummaryRefreshActor(override val serviceRegistryActor: ActorRef)
  extends LoggingFSM[SummaryRefreshState, SummaryRefreshData]
    with MetadataDatabaseAccess
    with MetadataServicesStore
    with CromwellInstrumentation {

  implicit val ec = context.dispatcher

  private val summaryMetricsGapsPath: NonEmptyList[String] = MetadataServiceActor.MetadataInstrumentationPrefix :+ "summarizer" :+ "gap"
  private val summaryMetricsProcessedPath: NonEmptyList[String] = MetadataServiceActor.MetadataInstrumentationPrefix :+ "summarizer" :+ "processed"

  val increasingGapPath = summaryMetricsGapsPath :+ "increasing"
  val decreasingGapPath = summaryMetricsGapsPath :+ "decreasing"

  val increasingProcessedPath = summaryMetricsProcessedPath :+ "increasing"
  val decreasingProcessedPath = summaryMetricsProcessedPath :+ "decreasing"

  private val instrumentationPrefix: Option[String] = InstrumentationPrefixes.ServicesPrefix

  startWith(WaitingForRequest, EmptySummaryRefreshData)

  when (WaitingForRequest) {
    case Event(SummarizeMetadata(limit, respondTo), data) =>

      val permittedSummaryStatusPointerUpdate: Option[Long] = data match {
        case EmptySummaryRefreshData => None
        case PreviousMaximumMetadataEntryId(value) => Option(value)
      }

      refreshWorkflowMetadataSummaries(limit, permittedSummaryStatusPointerUpdate) onComplete {
        case Success(summaryResult) =>
          sendGauge(increasingGapPath, summaryResult.increasingGap, instrumentationPrefix)
          sendGauge(decreasingGapPath, summaryResult.decreasingGap, instrumentationPrefix)

          count(increasingProcessedPath, summaryResult.rowsProcessedIncreasing, instrumentationPrefix)
          count(decreasingProcessedPath, summaryResult.rowsProcessedDecreasing, instrumentationPrefix)

          self ! MetadataSummaryComplete(PreviousMaximumMetadataEntryId(summaryResult.maximumKnownMetadataEntryId))
          respondTo ! MetadataSummarySuccess
        case Failure(t) =>
          log.error(t, "Failed to summarize metadata")
          self ! MetadataSummaryComplete(data)
          respondTo ! MetadataSummaryFailure(t)
      }
      goto(SummarizingMetadata)
  }

  when (SummarizingMetadata) {
    case Event(MetadataSummaryComplete(nextData: SummaryRefreshData), _) =>
      goto(WaitingForRequest) using nextData
  }

  whenUnhandled {
    case Event(wut, _) =>
      log.warning("Unrecognized or unexpected message while in state '{}': {}", stateName, wut)
      stay()
  }
}
