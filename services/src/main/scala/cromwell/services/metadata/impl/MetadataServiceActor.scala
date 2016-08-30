package cromwell.services.metadata.impl

import java.time.OffsetDateTime
import java.util.UUID

import akka.actor.{Actor, ActorLogging, Props}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
import cromwell.database.CromwellDatabase
import cromwell.services.metadata.MetadataService.{PutMetadataAction, ReadAction, RefreshSummary, ValidateWorkflowIdAndExecute}
import cromwell.services.metadata.impl.MetadataServiceActor._
import cromwell.services.metadata.impl.MetadataSummaryRefreshActor.{MetadataSummaryFailure, MetadataSummarySuccess, SummarizeMetadata}
import lenthall.config.ScalaConfig._

import scala.concurrent.duration.{Duration, FiniteDuration}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object MetadataServiceActor {

  val MetadataSummaryRefreshInterval =
    Duration(ConfigFactory.load().getStringOr("services.MetadataService.metadata-summary-refresh-interval", "2 seconds")).asInstanceOf[FiniteDuration]

  val MetadataSummaryTimestampMinimum =
    ConfigFactory.load().getStringOption("services.MetadataService.metadata-summary-timestamp-minimum") map OffsetDateTime.parse

  // A workflow will stay in the existence cache for this many runs of the workflow summary actor before being expired out.
  val CacheExpiryCount = 5

  def props(serviceConfig: Config, globalConfig: Config) = Props(MetadataServiceActor(serviceConfig, globalConfig))
}

case class MetadataServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with ActorLogging with MetadataDatabaseAccess with CromwellDatabase {

  val summaryActor = context.actorOf(MetadataSummaryRefreshActor.props(MetadataSummaryTimestampMinimum), "metadata-summary-actor")
  val readActor = context.actorOf(ReadMetadataActor.props(), "read-metadata-actor")
  val writeActor = context.actorOf(WriteMetadataActor.props(), "write-metadata-actor")
  implicit val ec = context.dispatcher

  self ! RefreshSummary

  private def scheduleSummary = context.system.scheduler.scheduleOnce(MetadataSummaryRefreshInterval, self, RefreshSummary)(context.dispatcher, self)

  private def validateWorkflowId(validation: ValidateWorkflowIdAndExecute): Unit = {
    val possibleWorkflowId = validation.possibleWorkflowId
    val callback = validation.validationCallback

    Try(UUID.fromString(possibleWorkflowId)) match {
      case Failure(t) => callback.onMalformed(possibleWorkflowId)
      case Success(uuid) =>
        workflowExistsWithId(possibleWorkflowId) onComplete {
          case Success(true) =>
            callback.onRecognized(WorkflowId(uuid))
          case Success(false) =>
            callback.onUnrecognized(possibleWorkflowId)
          case Failure(t) =>
            callback.onFailure(possibleWorkflowId, t)
        }
    }
  }

  def receive = {
    case action@PutMetadataAction(events) => writeActor forward action
    case v: ValidateWorkflowIdAndExecute => validateWorkflowId(v)
    case action: ReadAction => readActor forward action
    case RefreshSummary =>
      val sndr = sender()
      summaryActor ! SummarizeMetadata(sndr)
    case MetadataSummarySuccess => scheduleSummary
    case MetadataSummaryFailure(t) =>
      log.error(t, "Error summarizing metadata")
      scheduleSummary
  }
}
