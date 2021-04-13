package cromwell.services.metadata.impl.archiver

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.{WorkflowAborted, WorkflowFailed, WorkflowId, WorkflowSucceeded}
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus.Unarchived
import cromwell.services.metadata.WorkflowQueryKey._
import cromwell.services.metadata.impl.archiver.ArchiveMetadataSchedulerActor._
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}


class ArchiveMetadataSchedulerActor(archiveMetadataConfig: ArchiveMetadataConfig,
                                    override val serviceRegistryActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper with CromwellInstrumentation {

  implicit val ec: ExecutionContext = context.dispatcher

  // kick off archiving immediately
  self ! ArchiveNextWorkflowMessage

  override def receive: Receive = {
    case ArchiveNextWorkflowMessage => archiveNextWorkflow().onComplete({
      case Success(true) => self ! ArchiveNextWorkflowMessage
      case Success(false) => scheduleNextWorkflowToArchive()
      case Failure(error) =>
        log.error(error, s"Error while archiving, will retry.")
        scheduleNextWorkflowToArchive()
    })
    case ShutdownCommand => context.stop(self)  // TODO: cancel any streaming that might be happening
    case other => log.info(s"Programmer Error! The ArchiveMetadataSchedulerActor received unexpected message! ($sender sent $other})")
  }

  def archiveNextWorkflow(): Future[Boolean] = {
    for {
      maybeWorkflowId <- getNextWorkflowId()
      result <- maybeWorkflowId match {
        case Some(id) => for {
          _ <- streamMetadataToGcs(id)
          _ <- markWorkflowAsArchived(id)
        } yield true
        case None => Future.successful(false)
      }
    } yield result
  }

  def getNextWorkflowId(): Future[Option[WorkflowId]] = {
    ???
  }

  def streamMetadataToGcs(workflowId: WorkflowId): Future[Unit] = {
    ???
  }

  def markWorkflowAsArchived(workflowId: WorkflowId): Future[Unit] = {
    ???
  }

  def scheduleNextWorkflowToArchive(): Unit = {
    context.system.scheduler.scheduleOnce(archiveMetadataConfig.interval)(self ! ArchiveNextWorkflowMessage)
    ()
  }
}

object ArchiveMetadataSchedulerActor {
  case object ArchiveNextWorkflowMessage

  // TODO: Archive from oldest-first
  // TODO: Allow requirements like "End timestamp not within 1y (eg)"
  val queryParametersForWorkflowsToArchive: Seq[(String, String)] = Seq(
    IncludeSubworkflows.name -> "true",
    Status.name -> WorkflowSucceeded.toString,
    Status.name -> WorkflowFailed.toString,
    Status.name -> WorkflowAborted.toString,
    MetadataArchiveStatus.name -> Unarchived.toString,
    Page.name -> "1",
    PageSize.name -> "1"
  )

  def props(archiveMetadataConfig: ArchiveMetadataConfig, serviceRegistryActor: ActorRef): Props =
    Props(new ArchiveMetadataSchedulerActor(archiveMetadataConfig, serviceRegistryActor))
}
