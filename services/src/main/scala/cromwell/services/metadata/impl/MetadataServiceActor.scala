package cromwell.services.metadata.impl

import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.WorkflowId
import cromwell.services.SingletonServicesStore
import cromwell.services.metadata.MetadataService.{PutMetadataAction, ReadAction, RefreshSummary, ValidateWorkflowIdAndExecute}
import cromwell.services.metadata.impl.MetadataServiceActor._
import cromwell.services.metadata.impl.MetadataSummaryRefreshActor.{MetadataSummaryFailure, MetadataSummarySuccess, SummarizeMetadata}
import net.ceedubs.ficus.Ficus._
import scala.concurrent.duration.{Duration, FiniteDuration}
import scala.util.{Failure, Success, Try}

object MetadataServiceActor {

  val MetadataSummaryRefreshInterval: Option[FiniteDuration] = {
    val duration = Duration(ConfigFactory.load().as[Option[String]]("services.MetadataService.metadata-summary-refresh-interval").getOrElse("2 seconds"))
    if (duration.isFinite()) Option(duration.asInstanceOf[FiniteDuration]) else None
  }

  def props(serviceConfig: Config, globalConfig: Config) = Props(MetadataServiceActor(serviceConfig, globalConfig)).withDispatcher(ServiceDispatcher)
}

case class MetadataServiceActor(serviceConfig: Config, globalConfig: Config)
  extends Actor with ActorLogging with MetadataDatabaseAccess with SingletonServicesStore {

  private val summaryActor: Option[ActorRef] = buildSummaryActor

  val readActor = context.actorOf(ReadMetadataActor.props(), "read-metadata-actor")
  val writeActor = context.actorOf(WriteMetadataActor.props(), "write-metadata-actor")
  implicit val ec = context.dispatcher

  summaryActor foreach { _ => self ! RefreshSummary }

  private def scheduleSummary(): Unit = {
    MetadataSummaryRefreshInterval foreach { context.system.scheduler.scheduleOnce(_, self, RefreshSummary)(context.dispatcher, self) }
  }

  private def buildSummaryActor: Option[ActorRef] = {
    val actor = MetadataSummaryRefreshInterval map {
      _ => context.actorOf(MetadataSummaryRefreshActor.props(), "metadata-summary-actor")
    }
    val message = MetadataSummaryRefreshInterval match {
      case Some(interval) => s"Metadata summary refreshing every $interval."
      case None => "Metadata summary refresh is off."
    }
    log.info(message)
    actor
  }

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
    case RefreshSummary => summaryActor foreach { _ ! SummarizeMetadata(sender()) }
    case MetadataSummarySuccess => scheduleSummary()
    case MetadataSummaryFailure(t) =>
      log.error(t, "Error summarizing metadata")
      scheduleSummary()
  }
}
