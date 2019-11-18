package cromwell.engine.workflow.lifecycle.deletion

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.RootWorkflowId
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor._
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService.{GetRootAndSubworkflowOutputs, RootAndSubworkflowOutputsLookupFailure, RootAndSubworkflowOutputsLookupResponse, WorkflowOutputs, WorkflowOutputsFailure, WorkflowOutputsResponse}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.lang3.exception.ExceptionUtils

class DeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId,
                               serviceRegistryActor: ActorRef) extends LoggingFSM[DeleteWorkflowFilesActorState, DeleteWorkflowFilesActorStateData] {

//  override def receive: Receive = {
//    case StartWorkflowFilesDeletion =>
//      serviceRegistryActor ! GetRootAndSubworkflowOutputs(rootWorkflowId)
//    case o: RootAndSubworkflowOutputsLookupResponse =>
//      println(s"received list of outputs!! Outputs: $o")
//    case f: RootAndSubworkflowOutputsLookupFailure =>
//      println(s"Something went wrong. Error ${f.reason}")
//    case other =>
//      log.error(s"Programmer Error! The DeleteWorkflowFilesActor for root workflow ${rootWorkflowId.id.toString} " +
//      s"received unexpected message! ($sender sent $other})")
//  }

  startWith(Pending, NoData)

  when(Pending) {
    case Event(StartWorkflowFilesDeletion, NoData) => {
      serviceRegistryActor ! GetRootAndSubworkflowOutputs(rootWorkflowId)
      goto(FetchingAllOutputs) using NoData
    }
  }

  when(FetchingAllOutputs) {
    case Event(RootAndSubworkflowOutputsLookupResponse(_, metadataEvents), _) => {
      if (metadataEvents.nonEmpty) {
        serviceRegistryActor ! WorkflowOutputs(rootWorkflowId, convertResponseToJson = false)
        goto(FetchingFinalOutputs) using FetchingFinalOutputsData(metadataEvents)
      }
      else {
        log.info(s"Root workflow ${rootWorkflowId.id} does not have any outputs to delete.")
        stopSelf()
      }
    }
    case Event(RootAndSubworkflowOutputsLookupFailure(_, reason), _) => {
      log.error(s"Something went wrong while retrieving all outputs for root workflow ${rootWorkflowId.id} from database. " +
        s"Error: ${ExceptionUtils.getMessage(reason)}")
      stopSelf()
    }
  }

  when(FetchingFinalOutputs) {
    case Event(WorkflowOutputsResponse(_, metadataEvents), FetchingFinalOutputsData(allOutputs)) => {
      val intermediateOutputs = gatherIntermediateOutputs(allOutputs, metadataEvents)
      if (intermediateOutputs.nonEmpty) goto(DeletingIntermediateFiles) using DeletingIntermediateFilesData(intermediateOutputs)
      else {
        log.info(s"Root workflow ${rootWorkflowId.id} does not have any intermediate output files to delete.")
        stopSelf()
      }
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


  private def gatherIntermediateOutputs(allOutputsEvents: Seq[MetadataEvent], finalOutputsEvents: Seq[MetadataEvent]) = {
    def existsInFinalOutputs(metadataValue: String) = {
      finalOutputsEvents.map(p => p.value match {
        case Some(v) => v.value.equals(metadataValue)
        case None => false
      }).reduce(_ || _)
    }

//    if (finalOutputsEvents.nonEmpty) {
//      allOutputsEvents.filterNot(x => x.value match {
//        case Some(v) => existsInFinalOutputs(v.value)
//        case None => true
//      }).map(_.value.get.value).toSet
//    }
//    else allOutputsEvents.filter(x => x.value.nonEmpty).map(_.value.get.value).toSet

    if (finalOutputsEvents.nonEmpty) {
      allOutputsEvents.filterNot(x => x.value match {
        case Some(v) => existsInFinalOutputs(v.value)
        case None => true
      }).map(_.value.get.value).toSet
    }
    else allOutputsEvents.collect {
      case e if e.value.nonEmpty => e.value.get.value
    }.toSet
  }
}


object DeleteWorkflowFilesActor {

  // Commands
  sealed trait DeleteWorkflowFilesActorMessage
  object StartWorkflowFilesDeletion extends DeleteWorkflowFilesActorMessage

  sealed trait DeleteWorkflowFilesActorState
  object Pending extends DeleteWorkflowFilesActorState
  object FetchingAllOutputs extends DeleteWorkflowFilesActorState
  object FetchingFinalOutputs extends DeleteWorkflowFilesActorState
  object DeletingIntermediateFiles extends DeleteWorkflowFilesActorState

  sealed trait DeleteWorkflowFilesActorStateData
  object NoData extends DeleteWorkflowFilesActorStateData
  case class FetchingFinalOutputsData(allOutputs: Seq[MetadataEvent]) extends DeleteWorkflowFilesActorStateData
  case class DeletingIntermediateFilesData(intermediateFiles: Set[String]) extends DeleteWorkflowFilesActorStateData


  def props(rootWorkflowId: RootWorkflowId, serviceRegistryActor: ActorRef): Props = {
    Props(new DeleteWorkflowFilesActor(rootWorkflowId, serviceRegistryActor))
  }
}
