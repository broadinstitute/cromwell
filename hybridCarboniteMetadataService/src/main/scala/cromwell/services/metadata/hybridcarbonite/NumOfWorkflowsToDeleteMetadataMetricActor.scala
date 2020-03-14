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
import cromwell.services.metadata.hybridcarbonite.NumOfWorkflowsToDeleteMetadataMetricActor._
import cromwell.services.metadata.impl.{MetadataDatabaseAccess, MetadataServiceActor}

class NumOfWorkflowsToDeleteMetadataMetricActor(override val serviceRegistryActor: ActorRef)
  extends LoggingFSM[MetadataDeletionMetricsActorState, Unit]
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

  implicit val ec = context.dispatcher

  startWith(WaitingForMetricCalculationRequestOrMetricValue, ())

  when(WaitingForMetricCalculationRequestOrMetricValue) {
    case Event(NumOfWorkflowsToDeleteMetadataMetricValue(value), _) =>
      sendGauge(workflowsToDeleteMetadataMetricPath, value, InstrumentationPrefixes.ServicesPrefix)
      stay()
    case Event(CalculateNumOfWorkflowsToDeleteMetadataMetricValue(currentTimestampMinusDelay), _) =>
      countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
        MetadataArchiveStatus.toDatabaseValue(Archived),
        currentTimestampMinusDelay
      )
        .map(intVal => FinishedNumOfWorkflowsToDeleteMetadataMetricValueCalculation(intVal.toLong))
        .pipeTo(self)
      goto(MetricCalculationInProgress)
  }

  when(MetricCalculationInProgress) {
    case Event(CalculateNumOfWorkflowsToDeleteMetadataMetricValue, _) =>
      // metric actor is already busy calculating metric value, so we dismiss this request
      stay()
    case Event(FinishedNumOfWorkflowsToDeleteMetadataMetricValueCalculation(calculatedValue), _) =>
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

object NumOfWorkflowsToDeleteMetadataMetricActor {

  def props(serviceRegistryActor: ActorRef) = Props(new NumOfWorkflowsToDeleteMetadataMetricActor(serviceRegistryActor))

  private val workflowsToDeleteMetadataMetricPath: NonEmptyList[String] =
    MetadataServiceActor.MetadataInstrumentationPrefix :+ "delete" :+ "numOfWorkflowsToDeleteMetadata"

  // states
  sealed trait MetadataDeletionMetricsActorState
  case object WaitingForMetricCalculationRequestOrMetricValue extends MetadataDeletionMetricsActorState
  case object MetricCalculationInProgress extends MetadataDeletionMetricsActorState

  // messages
  case class CalculateNumOfWorkflowsToDeleteMetadataMetricValue(currentTimestampMinusDelay: OffsetDateTime)
  case class NumOfWorkflowsToDeleteMetadataMetricValue(value: Long)
  case class FinishedNumOfWorkflowsToDeleteMetadataMetricValueCalculation(calculatedValue: Long)

}
