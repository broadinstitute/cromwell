package cromwell.services.metadata

import java.time.OffsetDateTime
import java.util.UUID

import akka.actor.{Actor, ActorLogging}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
import cromwell.engine.db.DataAccess
import cromwell.services.MetadataServiceActor._
import MetadataSummaryRefreshActor.{MetadataSummaryFailure, MetadataSummarySuccess, SummarizeMetadata}
import cromwell.services.metadata.EngineMetadataServiceActor._
import lenthall.config.ScalaConfig._

import scala.concurrent.duration.{Duration, FiniteDuration}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}


object EngineMetadataServiceActor {

  val MetadataSummaryRefreshInterval =
    Duration(ConfigFactory.load().getStringOr("services.MetadataService.metadata-summary-refresh-interval", "2 seconds")).asInstanceOf[FiniteDuration]

  val MetadataSummaryTimestampMinimum =
    ConfigFactory.load().getStringOption("services.MetadataService.metadata-summary-timestamp-minimum") map OffsetDateTime.parse

  // A workflow will stay in the existence cache for this many runs of the workflow summary actor before being expired out.
  val CacheExpiryCount = 5
}

// TODO: PBE: Will not be MetadataServiceActor until circular dependencies fixed.
case class EngineMetadataServiceActor(serviceConfig: Config, globalConfig: Config) extends Actor with ActorLogging {

  val dataAccess = DataAccess.globalDataAccess
  val summaryActor = context.actorOf(MetadataSummaryRefreshActor.props(MetadataSummaryTimestampMinimum), "metadata-summary-actor")
  val readActor = context.actorOf(ReadMetadataActor.props(), "read-metadata-actor")
  val writeActor = context.actorOf(WriteMetadataActor.props(), "write-metadata-actor")
  implicit val ec = context.dispatcher

  self ! RefreshSummary

  private def scheduleSummary = context.system.scheduler.scheduleOnce(MetadataSummaryRefreshInterval, self, RefreshSummary)(context.dispatcher)

  private def validateWorkflowId(validation: ValidateWorkflowIdAndExecute): Unit = {
    val possibleWorkflowId = validation.possibleWorkflowId
    val requestContext = validation.requestContext
    val callback = validation.validationCallback

    Try(UUID.fromString(possibleWorkflowId)) match {
      case Failure(t) => callback.onMalformed(possibleWorkflowId)(requestContext)
      case Success(uuid) =>
        dataAccess.workflowExistsWithId(possibleWorkflowId) onComplete {
          case Success(true) =>
            callback.onRecognized(WorkflowId(uuid))(requestContext)
          case Success(false) =>
            callback.onUnrecognized(possibleWorkflowId)(requestContext)
          case Failure(t) =>
            callback.onFailure(possibleWorkflowId, t)(requestContext)
        }
    }
  }

  def receive = {
    case action@PutMetadataAction(events) => writeActor forward action
    case v: ValidateWorkflowIdAndExecute => validateWorkflowId(v)
    case action: ReadAction => readActor forward action
    case RefreshSummary => summaryActor ! SummarizeMetadata
    case MetadataSummarySuccess => scheduleSummary
    case MetadataSummaryFailure(t) =>
      log.error(t, "Error summarizing metadata")
      scheduleSummary
  }
}
