package cromwell.services.metadata.impl.archiver

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.typesafe.config.Config
import common.util.TimeUtil.EnhancedOffsetDateTime
import cromwell.services.MetadataServicesStore
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metadata.impl.MetadataDatabaseAccess

import scala.concurrent.duration.FiniteDuration

/*
  This class would do same things cromwell.services.metadata.hybridcarbonite.DeleteMetadataActor does
 */

class DeleteMetadataActor(metadataDeletionConfig: Config,
                          override val serviceRegistryActor: ActorRef)
  extends Actor
  with ActorLogging
  with MetadataDatabaseAccess
  with MetadataServicesStore
  with CromwellInstrumentation {

  override def receive: Receive = ???
}

object DeleteMetadataActor {
  import cromwell.services.metadata.WorkflowQueryKey._
  import cromwell.services.metadata.MetadataArchiveStatus.Archived
  import cromwell.core.{WorkflowAborted, WorkflowFailed, WorkflowSucceeded}

  import java.time.OffsetDateTime

  def queryParametersForWorkflowsToDelete(currentTime: OffsetDateTime, deleteDelay: FiniteDuration): Seq[(String, String)] = Seq(
    IncludeSubworkflows.name -> "true",
    Status.name -> WorkflowSucceeded.toString,
    Status.name -> WorkflowFailed.toString,
    Status.name -> WorkflowAborted.toString,
    MetadataArchiveStatus.name -> Archived.toString,
    Page.name -> "1",
    PageSize.name -> "1",
    NewestFirst.name -> "false", // oldest first for deleting
    EndDate.name -> currentTime.minusNanos(deleteDelay.toNanos).toUtcMilliString
  )
}
