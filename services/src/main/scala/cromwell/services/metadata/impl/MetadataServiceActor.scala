package cromwell.services.metadata.impl


import akka.actor.SupervisorStrategy.{Decider, Directive, Escalate, Resume}
import akka.actor.{Actor, ActorContext, ActorInitializationException, ActorLogging, ActorRef, Cancellable, OneForOneStrategy, Props}
import akka.routing.Listen
import cats.data.NonEmptyList
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.{LoadConfig, WorkflowId}
import cromwell.services.MetadataServicesStore
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.impl.MetadataServiceActor._
import cromwell.services.metadata.impl.MetadataSummaryRefreshActor.{MetadataSummaryFailure, MetadataSummarySuccess, SummarizeMetadata}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

object MetadataServiceActor {
  val MetadataInstrumentationPrefix = NonEmptyList.of("metadata")

  val MetadataSummaryRefreshInterval: Option[FiniteDuration] = {
    val duration = Duration(ConfigFactory.load().as[Option[String]]("services.MetadataService.config.metadata-summary-refresh-interval").getOrElse("2 seconds"))
    if (duration.isFinite()) Option(duration.asInstanceOf[FiniteDuration]) else None
  }

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(MetadataServiceActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
}

final case class MetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
  extends Actor with ActorLogging with MetadataDatabaseAccess with MetadataServicesStore with GracefulShutdownHelper {

  private val decider: Decider = {
    case _: ActorInitializationException => Escalate
    case _ => Resume
  }

  override val supervisorStrategy = new OneForOneStrategy()(decider) {
    override def logFailure(context: ActorContext, child: ActorRef, cause: Throwable, decision: Directive) = {
      val childName = if (child == readActor) "Read" else "Write"
      log.error(s"The $childName Metadata Actor died unexpectedly, metadata events might have been lost. Restarting it...", cause)
    }
  }

  private val summaryActor: Option[ActorRef] = buildSummaryActor

  val readActor = context.actorOf(ReadMetadataActor.props(), "read-metadata-actor")

  val dbFlushRate = serviceConfig.as[Option[FiniteDuration]]("db-flush-rate").getOrElse(5 seconds)
  val dbBatchSize = serviceConfig.as[Option[Int]]("db-batch-size").getOrElse(200)
  val writeActor = context.actorOf(WriteMetadataActor.props(dbBatchSize, dbFlushRate, serviceRegistryActor, LoadConfig.MetadataWriteThreshold), "WriteMetadataActor")
  implicit val ec = context.dispatcher
  private var summaryRefreshCancellable: Option[Cancellable] = None

  summaryActor foreach { _ => self ! RefreshSummary }

  private def scheduleSummary(): Unit = {
    MetadataSummaryRefreshInterval foreach { interval =>
      summaryRefreshCancellable = Option(context.system.scheduler.scheduleOnce(interval, self, RefreshSummary)(context.dispatcher, self))
    }
  }

  override def postStop(): Unit = {
    summaryRefreshCancellable foreach { _.cancel() }
    super.postStop()
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

  private def validateWorkflowIdInMetadata(possibleWorkflowId: WorkflowId, sender: ActorRef): Unit = {
    workflowWithIdExistsInMetadata(possibleWorkflowId.toString) onComplete {
      case Success(true) => sender ! RecognizedWorkflowId
      case Success(false) => sender ! UnrecognizedWorkflowId
      case Failure(e) => sender ! FailedToCheckWorkflowId(new RuntimeException(s"Failed lookup attempt for workflow ID $possibleWorkflowId", e))
    }
  }

  private def validateWorkflowIdInMetadataSummaries(possibleWorkflowId: WorkflowId, sender: ActorRef): Unit = {
    workflowWithIdExistsInMetadataSummaries(possibleWorkflowId.toString) onComplete {
      case Success(true) => sender ! RecognizedWorkflowId
      case Success(false) => sender ! UnrecognizedWorkflowId
      case Failure(e) => sender ! FailedToCheckWorkflowId(new RuntimeException(s"Failed lookup attempt for workflow ID $possibleWorkflowId", e))
    }
  }

  def receive = {
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(writeActor))
    case action: PutMetadataAction => writeActor forward action
    case action: PutMetadataActionAndRespond => writeActor forward action
    // Assume that listen messages are directed to the write metadata actor
    case listen: Listen => writeActor forward listen
    case v: ValidateWorkflowIdInMetadata => validateWorkflowIdInMetadata(v.possibleWorkflowId, sender())
    case v: ValidateWorkflowIdInMetadataSummaries => validateWorkflowIdInMetadataSummaries(v.possibleWorkflowId, sender())
    case action: ReadAction => readActor forward action
    case RefreshSummary => summaryActor foreach { _ ! SummarizeMetadata(sender()) }
    case MetadataSummarySuccess => scheduleSummary()
    case MetadataSummaryFailure(t) =>
      log.error(t, "Error summarizing metadata")
      scheduleSummary()
  }
}
