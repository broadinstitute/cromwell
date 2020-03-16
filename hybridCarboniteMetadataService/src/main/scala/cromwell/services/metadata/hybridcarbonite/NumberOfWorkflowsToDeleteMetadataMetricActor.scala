package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.actor.{ActorRef, LoggingFSM, Props, Status}
import akka.pattern.pipe
import cats.data.NonEmptyList
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus.Archived
import cromwell.services.metadata.hybridcarbonite.NumberOfWorkflowsToDeleteMetadataMetricActor._
import cromwell.services.metadata.impl.{MetadataDatabaseAccess, MetadataServiceActor}

class NumberOfWorkflowsToDeleteMetadataMetricActor(override val serviceRegistryActor: ActorRef)
  extends LoggingFSM[MetadataDeletionMetricsActorState, Unit]
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

  implicit val ec = context.dispatcher

  startWith(WaitingForMetricCalculationRequestOrMetricValue, ())

  when(WaitingForMetricCalculationRequestOrMetricValue) {
    case Event(NumberOfWorkflowsToDeleteMetadataMetricValue(value), _) =>
      sendGauge(workflowsToDeleteMetadataMetricPath, value, InstrumentationPrefixes.ServicesPrefix)
      stay()
    case Event(CalculateNumberOfWorkflowsToDeleteMetadataMetricValue(currentTimestampMinusDelay), _) =>
      countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
        MetadataArchiveStatus.toDatabaseValue(Archived),
        currentTimestampMinusDelay
      )
        .map(intVal => FinishedNumberOfWorkflowsToDeleteMetadataMetricValueCalculation(intVal.toLong))
        .pipeTo(self)
      goto(MetricCalculationInProgress)
  }

  when(MetricCalculationInProgress) {
    case Event(CalculateNumberOfWorkflowsToDeleteMetadataMetricValue, _) =>
      // metric actor is already busy calculating metric value, so we dismiss this request
      stay()
    case Event(NumberOfWorkflowsToDeleteMetadataMetricValue(value), _) =>
      // pass through the pre-calculated metric value and stay, continuing to wait for DB query to finish
      sendGauge(workflowsToDeleteMetadataMetricPath, value, InstrumentationPrefixes.ServicesPrefix)
      stay()
    case Event(FinishedNumberOfWorkflowsToDeleteMetadataMetricValueCalculation(calculatedValue), _) =>
      // populate metric and return to the default state
      sendGauge(workflowsToDeleteMetadataMetricPath, calculatedValue, InstrumentationPrefixes.ServicesPrefix)
      goto(WaitingForMetricCalculationRequestOrMetricValue)
    case Event(Status.Failure(ex), _) =>
      log.error(s"Cannot send numOfworkflowsToDeleteMetadata gauge metric value: unable to query the value from the database.", ex)
      goto(WaitingForMetricCalculationRequestOrMetricValue)
  }

  whenUnhandled {
    case Event(unexpected, _) =>
      log.warning(s"Programmer error: this actor should not receive message $unexpected while in state $stateName")
      stay()
  }
}

object NumberOfWorkflowsToDeleteMetadataMetricActor {

  def props(serviceRegistryActor: ActorRef) = Props(new NumberOfWorkflowsToDeleteMetadataMetricActor(serviceRegistryActor))

  private val workflowsToDeleteMetadataMetricPath: NonEmptyList[String] =
    MetadataServiceActor.MetadataInstrumentationPrefix :+ "delete" :+ "numOfWorkflowsToDeleteMetadata"

  // states
  sealed trait MetadataDeletionMetricsActorState
  case object WaitingForMetricCalculationRequestOrMetricValue extends MetadataDeletionMetricsActorState
  case object MetricCalculationInProgress extends MetadataDeletionMetricsActorState

  // messages
  case class CalculateNumberOfWorkflowsToDeleteMetadataMetricValue(currentTimestampMinusDelay: OffsetDateTime)
  case class NumberOfWorkflowsToDeleteMetadataMetricValue(value: Long)
  case class FinishedNumberOfWorkflowsToDeleteMetadataMetricValueCalculation(calculatedValue: Long)

}
