package cromwell.engine.workflow.lifecycle.deletion

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.RootWorkflowId
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder, PathFactory}
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor._
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService.{GetRootAndSubworkflowOutputs, RootAndSubworkflowOutputsLookupFailure, RootAndSubworkflowOutputsLookupResponse, WorkflowOutputs, WorkflowOutputsFailure, WorkflowOutputsResponse}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.lang3.exception.ExceptionUtils

import scala.util.Try

class DeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId,
                               serviceRegistryActor: ActorRef,
                               pathBuilders: List[PathBuilder]) extends LoggingFSM[DeleteWorkflowFilesActorState, DeleteWorkflowFilesActorStateData] {

  /*
  building a path for a string output such as 'Hello World' satisfies the conditions of a
  being a relative file path for a DefaultPath which we don't want, so remove DefaultPathBuilder
  (this is because we only want cloud backend paths like gs://, s3://, etc.)
   */
  val pathBuildersWithoutDefault = pathBuilders.filterNot(p => p == DefaultPathBuilder)

  startWith(Pending, NoData)

  when(Pending) {
    case Event(StartWorkflowFilesDeletion, NoData) =>
      serviceRegistryActor ! GetRootAndSubworkflowOutputs(rootWorkflowId)
      goto(FetchingAllOutputs) using NoData
  }

  when(FetchingAllOutputs) {
    case Event(RootAndSubworkflowOutputsLookupResponse(_, metadataEvents), _) =>
      if (metadataEvents.nonEmpty) {
        serviceRegistryActor ! WorkflowOutputs(rootWorkflowId, convertResponseToJson = false)
        goto(FetchingFinalOutputs) using FetchingFinalOutputsData(metadataEvents)
      }
      else {
        log.info(s"Root workflow ${rootWorkflowId.id} does not have any outputs to delete.")
        stopSelf()
      }
    case Event(RootAndSubworkflowOutputsLookupFailure(_, reason), _) =>
      log.error(s"Something went wrong while retrieving all outputs for root workflow ${rootWorkflowId.id} from database. " +
        s"Error: ${ExceptionUtils.getMessage(reason)}")
      stopSelf()
  }

  when(FetchingFinalOutputs) {
    case Event(WorkflowOutputsResponse(_, metadataEvents), FetchingFinalOutputsData(allOutputs)) =>
      val intermediateOutputs = gatherIntermediateOutputFiles(allOutputs, metadataEvents)
      if (intermediateOutputs.nonEmpty) goto(DeletingIntermediateFiles) using DeletingIntermediateFilesData(intermediateOutputs)
      else {
        log.info(s"Root workflow ${rootWorkflowId.id} does not have any intermediate output files to delete.")
        stopSelf()
      }
    case Event(WorkflowOutputsFailure(_, reason), _) =>
      log.error(s"Something went wrong while retrieving final outputs for root workflow ${rootWorkflowId.id} from database. " +
        s"Error: ${ExceptionUtils.getMessage(reason)}")
      stopSelf()
  }

  //TODO: To be done as part of WA-41 (https://broadworkbench.atlassian.net/browse/WA-41)
  when(DeletingIntermediateFiles) {
    case Event(_, _) => stopSelf()
  }

  whenUnhandled {
    case Event(ShutdownCommand, _) => stopSelf()
    case other =>
      log.error(s"Programmer Error: Unexpected message to ${getClass.getSimpleName} ${self.path.name} in state $stateName with $stateData: $other")
      stay()
  }


  private def stopSelf() = {
    context stop self
    stay()
  }

  def gatherIntermediateOutputFiles(allOutputsEvents: Seq[MetadataEvent], finalOutputsEvents: Seq[MetadataEvent]): Set[Path] = {
    val isFinalOutputsNonEmpty = finalOutputsEvents.nonEmpty

    def existsInFinalOutputs(metadataValue: String): Boolean = {
      finalOutputsEvents.map(p => p.value match {
        case Some(v) => v.value.equals(metadataValue)
        case None => false
      }).reduce(_ || _)
    }

    def getFilePath(metadataValue: String): Option[Path] = {
      //if workflow has final outputs check if this output value does not exist in final outputs
      (isFinalOutputsNonEmpty, existsInFinalOutputs(metadataValue)) match {
        case (true, true) => None
        case (true, false) | (false, _) =>
          // eliminate outputs which are not files
          Try(PathFactory.buildPath(metadataValue, pathBuildersWithoutDefault)).toOption
      }
    }

    allOutputsEvents.collect {
      case o if o.value.nonEmpty => getFilePath(o.value.get.value)
    }.flatten.toSet
  }
}


object DeleteWorkflowFilesActor {

  // Commands
  sealed trait DeleteWorkflowFilesActorMessage
  object StartWorkflowFilesDeletion extends DeleteWorkflowFilesActorMessage

  // Actor states
  sealed trait DeleteWorkflowFilesActorState
  object Pending extends DeleteWorkflowFilesActorState
  object FetchingAllOutputs extends DeleteWorkflowFilesActorState
  object FetchingFinalOutputs extends DeleteWorkflowFilesActorState
  object DeletingIntermediateFiles extends DeleteWorkflowFilesActorState

  // State data
  sealed trait DeleteWorkflowFilesActorStateData
  object NoData extends DeleteWorkflowFilesActorStateData
  case class FetchingFinalOutputsData(allOutputs: Seq[MetadataEvent]) extends DeleteWorkflowFilesActorStateData
  case class DeletingIntermediateFilesData(intermediateFiles: Set[Path]) extends DeleteWorkflowFilesActorStateData


  def props(rootWorkflowId: RootWorkflowId, serviceRegistryActor: ActorRef, pathBuilders: List[PathBuilder]): Props = {
    Props(new DeleteWorkflowFilesActor(rootWorkflowId, serviceRegistryActor, pathBuilders))
  }
}
