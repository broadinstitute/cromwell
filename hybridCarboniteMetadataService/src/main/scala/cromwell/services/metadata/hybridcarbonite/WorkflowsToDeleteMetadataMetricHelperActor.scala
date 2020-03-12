package cromwell.services.metadata.hybridcarbonite

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus.Archived
import cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor.MetricActorFreed
import cromwell.services.metadata.hybridcarbonite.WorkflowsToDeleteMetadataMetricHelperActor._
import cromwell.services.metadata.impl.{MetadataDatabaseAccess, MetadataServiceActor}

import scala.util.{Failure, Success}

class WorkflowsToDeleteMetadataMetricHelperActor(metadataDeletionActor: ActorRef, override val serviceRegistryActor: ActorRef)
  extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

  implicit val ec = context.dispatcher

  override def postRestart(reason: Throwable): Unit = {
    super.postRestart(reason)
    metadataDeletionActor ! MetricActorFreed
  }

  override def receive: Receive = {
    case SendPrecalculatedMetric(workflowsToDeleteMetadata) =>
      sendGauge(workflowsToDeleteMetadataMetricPath, workflowsToDeleteMetadata, InstrumentationPrefixes.ServicesPrefix)
      metadataDeletionActor ! MetricActorFreed
    case CalculateMetricAndSend(currentTimestampMinusDelay) =>
      countRootWorkflowSummaryEntriesByArchiveStatusAndOlderThanTimestamp(
        MetadataArchiveStatus.toDatabaseValue(Archived),
        currentTimestampMinusDelay
      ) onComplete {
        case Success(numToDelete) =>
          sendGauge(workflowsToDeleteMetadataMetricPath, numToDelete.toLong, InstrumentationPrefixes.ServicesPrefix)
          metadataDeletionActor ! MetricActorFreed
        case Failure(ex) =>
          log.error(s"Cannot send workflowsToDeleteMetadata gauge metric value: unable to query number of " +
            s"workflows-to-delete-metadata from the database.", ex)
          metadataDeletionActor ! MetricActorFreed
      }
  }

}

object WorkflowsToDeleteMetadataMetricHelperActor {

  def props(metadataDeletionActor: ActorRef, serviceRegistryActor: ActorRef) =
    Props(new WorkflowsToDeleteMetadataMetricHelperActor(metadataDeletionActor, serviceRegistryActor))

  case class CalculateMetricAndSend(currentTimestampMinusDelay: OffsetDateTime)
  case class SendPrecalculatedMetric(value: Long)

  private val workflowsToDeleteMetadataMetricPath: NonEmptyList[String] =
    MetadataServiceActor.MetadataInstrumentationPrefix :+ "delete" :+ "workflowsToDeleteMetadata"

}
