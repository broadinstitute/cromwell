package cromwell.engine.workflow.lifecycle.deletion

import akka.actor.{LoggingFSM, Props}
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder, PathFactory}
import cromwell.core.{CallOutputs, RootWorkflowId}
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wom.values.{WomSingleFile, WomValue}

import scala.util.Try

class DeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId,
                               workflowFinalOutputs: CallOutputs,
                               workflowAllOutputs: CallOutputs,
                               pathBuilders: List[PathBuilder]) extends LoggingFSM[DeleteWorkflowFilesActorState, DeleteWorkflowFilesActorStateData] {

  /*
  building a path for a string output such as 'Hello World' satisfies the conditions of a
  being a relative file path for a DefaultPath which we don't want, so remove DefaultPathBuilder
  (this is because we only want cloud backend paths like gs://, s3://, etc.)
   */
  val pathBuildersWithoutDefault = pathBuilders.filterNot(_ == DefaultPathBuilder)

  startWith(Pending, NoData)

  when(Pending) {
    case Event(StartWorkflowFilesDeletion, NoData) =>
      val intermediateOutputs = gatherIntermediateOutputFiles(workflowAllOutputs.outputs, workflowFinalOutputs.outputs)
      if (intermediateOutputs.nonEmpty) goto(DeletingIntermediateFiles) using DeletingIntermediateFilesData(intermediateOutputs)
      else {
        log.info(s"Root workflow ${rootWorkflowId.id} does not have any intermediate output files to delete.")
        stopSelf()
      }
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

  def gatherIntermediateOutputFiles(allOutputs: Map[OutputPort, WomValue], finalOutputs: Map[OutputPort, WomValue]): Set[Path] = {
    val isFinalOutputsNonEmpty = finalOutputs.nonEmpty

    def existsInFinalOutputs(value: String): Option[Boolean] = {
      finalOutputs.map(p => value.equals(p._2.valueString)).reduceOption(_ || _)
    }

    def getFilePath(value: String): Option[Path] = {
      //if workflow has final outputs check if this output value does not exist in final outputs
      (isFinalOutputsNonEmpty, existsInFinalOutputs(value)) match {
        case (true, Some(true)) => None
        case (true, Some(false)) | (true, None) | (false, _) =>
          // eliminate outputs which are not files
          Try(PathFactory.buildPath(value, pathBuildersWithoutDefault)).toOption
      }
    }

    def checkOutputIsFile(port: OutputPort, value: WomValue): Boolean = {
      (port, value) match {
        case (_: GraphNodeOutputPort, _: WomSingleFile) => true
        case _ => false
      }
    }


    allOutputs.collect {
      case o if checkOutputIsFile(o._1, o._2) => getFilePath(o._2.valueString)
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
  object DeletingIntermediateFiles extends DeleteWorkflowFilesActorState

  // State data
  sealed trait DeleteWorkflowFilesActorStateData
  object NoData extends DeleteWorkflowFilesActorStateData
  case class DeletingIntermediateFilesData(intermediateFiles: Set[Path]) extends DeleteWorkflowFilesActorStateData


  def props(rootWorkflowId: RootWorkflowId,
            workflowFinalOutputs: CallOutputs,
            workflowAllOutputs: CallOutputs,
            pathBuilders: List[PathBuilder]): Props = {
    Props(new DeleteWorkflowFilesActor(rootWorkflowId, workflowFinalOutputs, workflowAllOutputs, pathBuilders))
  }
}
