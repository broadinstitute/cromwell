package cromwell.services.metadata.impl.deleter

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.util.Timeout
import common.util.StringUtil.EnhancedToStringable
import cromwell.core.WorkflowId
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataArchiveStatus.{Archived, ArchivedAndPurged}
import cromwell.services.metadata.impl.MetadataDatabaseAccess
import cromwell.services.metadata.impl.deleter.DeleteMetadataActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

class DeleteMetadataActor(deleteMetadataConfig: DeleteMetadataConfig,
                          override val serviceRegistryActor: ActorRef)
  extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

  implicit val ec: ExecutionContext = context.dispatcher
  implicit val askTimeout: Timeout = new Timeout(60.seconds)

  // Send an initial delete message to get things started:
  self ! DeleteNextWorkflowMessage

  override def receive: Receive = {
    case DeleteNextWorkflowMessage => deleteNextWorkflow().onComplete({
      case Success(true) => self ! DeleteNextWorkflowMessage
      case Success(false) => scheduleNextDeleteAttemptAfterInterval()
      case Failure(error) =>
        log.error(error, s"Error while deleting, will wait ${deleteMetadataConfig.nonSuccessInterval} then try again.")
        scheduleNextDeleteAttemptAfterInterval()
    })
    case ShutdownCommand => context.stop(self)  // TODO: cancel any deletion action that might be happening?
    case other => log.info(s"Programmer Error! The ArchiveMetadataSchedulerActor received unexpected message! ($sender sent ${other.toPrettyElidedString(1000)}})")
  }

  def deleteNextWorkflow(): Future[Boolean] = for {
    maybeWorkflowId <- lookupNextWorkflowToDelete()
    result <- maybeWorkflowId match {
      case Some(id) =>
        log.info(s"Workflow $id identified for metadata deletion")
        for {
          rowsDeleted <- deleteAllMetadataEntriesForWorkflowAndUpdateArchiveStatus(id, MetadataArchiveStatus.toDatabaseValue(ArchivedAndPurged))
          _ = log.info(s"Deleted $rowsDeleted metadata rows for $id")
        } yield true
      case None => Future.successful(false)
    }
  } yield result

  def lookupNextWorkflowToDelete(): Future[Option[WorkflowId]] = {
    val currentTimestampMinusDelay = OffsetDateTime.now().minusSeconds(deleteMetadataConfig.delayAfterWorkflowCompletion.toSeconds)

    for {
      queryResult <- queryWorkflowIdsByArchiveStatusAndOlderThanTimestamp(
        MetadataArchiveStatus.toDatabaseValue(Archived),
        currentTimestampMinusDelay,
        batchSize = 1
      )
      headOption = queryResult.headOption
      workflowId = headOption.map(WorkflowId.fromString)
    } yield workflowId
  }

  def scheduleNextDeleteAttemptAfterInterval(): Unit = {
    context.system.scheduler.scheduleOnce(deleteMetadataConfig.nonSuccessInterval)(self ! DeleteNextWorkflowMessage)
    ()
  }

}

object DeleteMetadataActor {
  case object DeleteNextWorkflowMessage

  def props(config: DeleteMetadataConfig, serviceRegistryActor: ActorRef): Props = Props(new DeleteMetadataActor(config, serviceRegistryActor))
}
