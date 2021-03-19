package cromwell.services.metadata.impl.archiver

import akka.actor.{Actor, ActorLogging, Props}
import cromwell.core.WorkflowId
import cromwell.services.metadata.impl.archiver.StreamMetadataToGcsActor.{ArchiveMetadataForWorkflow, MetadataArchiveResult}
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.lang3.NotImplementedException

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.{Failure, Try}

/*
  This actor will archive the metadata by streaming it to GCS. Almost similar to what
  cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor does
 */
class StreamMetadataToGcsActor extends Actor with ActorLogging with GracefulShutdownHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case ArchiveMetadataForWorkflow(_) =>
      val s = sender()
      context.system.scheduler.scheduleOnce(2.minutes) {
        s ! MetadataArchiveResult(Failure(new NotImplementedException("Archive process not implemented yet!")))
      }
    case ShutdownCommand => context.stop(self) // TODO: Anything more graceful?
    case other => log.info(s"Programmer Error! The ArchiveMetadataSchedulerActor received unexpected message! ($sender sent $other})")
  }
}

object StreamMetadataToGcsActor {
  def props(): Props = Props(new StreamMetadataToGcsActor)

  sealed trait StreamMetadataToGcsActorMessage
  final case class ArchiveMetadataForWorkflow(workflow: WorkflowId) extends StreamMetadataToGcsActorMessage
  final case class MetadataArchiveResult(freezingResult: Try[Unit]) extends StreamMetadataToGcsActorMessage

}