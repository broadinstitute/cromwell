package cromwell.services.metadata.impl.archiver

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

/*
  This class would do what CarboniteWorkerActor is doing. It would schedule workflows to archive at periodic intervals
  and call StreamMetadataToGcsActor, who shall actually stream metadata to GCS
 */
class ArchiveMetadataSchedulerActor(archiveMetadataConfig: ArchiveMetadataConfig,
                                    override val serviceRegistryActor: ActorRef) extends Actor with ActorLogging with GracefulShutdownHelper with CromwellInstrumentation {

  // val streamMetadataActor = ??? // instance of StreamMetadataToGcsActor

  override def receive: Receive = {
    case ShutdownCommand => context.stop(self) // TODO: handle this better when we are actually archiving
    case other => log.info(s"Programmer Error! The ArchiveMetadataSchedulerActor received unexpected message! ($sender sent $other})")
  }
}


object ArchiveMetadataSchedulerActor {

  def props(archiveMetadataConfig: ArchiveMetadataConfig, serviceRegistryActor: ActorRef) =
    Props(new ArchiveMetadataSchedulerActor(archiveMetadataConfig, serviceRegistryActor))
}
