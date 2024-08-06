package cromwell.engine.workflow.lifecycle.finalization

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.event.LoggingReceive
import cromwell.backend.BackendLifecycleActor.BackendWorkflowLifecycleActorResponse
import cromwell.backend.BackendWorkflowFinalizationActor.{
  FinalizationFailed,
  FinalizationResponse,
  FinalizationSuccess,
  Finalize
}
import cromwell.backend.AllBackendInitializationData
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.WorkflowOptions.{Copy, Move}
import cromwell.core._
import cromwell.core.io.AsyncIoActorClient
import cromwell.core.path.Path
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.OutputsLocationHelper
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object CopyWorkflowOutputsActor {
  def props(workflowId: WorkflowId,
            ioActor: ActorRef,
            workflowDescriptor: EngineWorkflowDescriptor,
            workflowOutputs: CallOutputs,
            initializationData: AllBackendInitializationData
  ) = Props(
    new CopyWorkflowOutputsActor(workflowId, ioActor, workflowDescriptor, workflowOutputs, initializationData)
  ).withDispatcher(IoDispatcher)
}

class CopyWorkflowOutputsActor(workflowId: WorkflowId,
                               override val ioActor: ActorRef,
                               val workflowDescriptor: EngineWorkflowDescriptor,
                               workflowOutputs: CallOutputs,
                               initializationData: AllBackendInitializationData
) extends Actor
    with ActorLogging
    with AsyncIoActorClient
    with OutputsLocationHelper {
  override lazy val ioCommandBuilder = GcsBatchCommandBuilder
  implicit val ec = context.dispatcher

  override def receive = LoggingReceive { case Finalize =>
    performActionThenRespond(afterAll()(context.dispatcher), FinalizationFailed)(context.dispatcher)
  }

  private def performActionThenRespond(operation: => Future[BackendWorkflowLifecycleActorResponse],
                                       onFailure: (Throwable) => BackendWorkflowLifecycleActorResponse
  )(implicit ec: ExecutionContext) = {
    val respondTo: ActorRef = sender()
    operation onComplete {
      case Success(r) => respondTo ! r
      case Failure(t) => respondTo ! onFailure(t)
    }
  }

  private def markDuplicates(outputFilePaths: Map[Path, Path]) = {
    // Check if there are duplicated destination paths and throw an exception if that is the case.
    // This creates a map of destinations and source paths which point to them in cases where there are multiple
    // source paths that point to the same destination.
    val duplicatedDestPaths: Map[Path, List[Path]] =
      outputFilePaths.groupBy { case (_, destPath) => destPath }.collect {
        case (destPath, list) if list.size > 1 => destPath -> list.toList.map { case (source, _) => source }
      }
    if (duplicatedDestPaths.nonEmpty) {
      val formattedCollidingCopyOptions = duplicatedDestPaths.toList
        .sortBy { case (dest, _) => dest.pathAsString } // Sort by destination path
        // Make a '/my/src -> /my/dest' copy tape string for each source and destination. Use flat map to get a single list
        // srcList is also sorted to get a deterministic output order. This is necessary for making sure the tests
        // for the error always succeed.
        .flatMap { case (dest, srcList) => srcList.sortBy(_.pathAsString).map(_.pathAsString + s" -> $dest") }
      throw new IllegalStateException(
        "Cannot copy output files to given final_workflow_outputs_dir" +
          s" as multiple files will be copied to the same path: \n${formattedCollidingCopyOptions.mkString("\n")}"
      )
    }
  }

  private def copyWorkflowOutputs(outputsDir: String): Future[Seq[Unit]] = {
    val outputFilePaths =
      outputFilePathMapping(outputsDir, workflowDescriptor, initializationData, workflowOutputs.outputs.values.toSeq)

    markDuplicates(outputFilePaths)

    val copies = outputFilePaths.toList map { case (srcPath, dstPath) =>
      asyncIo.copyAsync(srcPath, dstPath)
    }

    Future.sequence(copies)
  }

  private def moveWorkflowOutputs(outputsDir: String): Future[Seq[Unit]] = {
    val outputFilePaths =
      outputFilePathMapping(outputsDir, workflowDescriptor, initializationData, workflowOutputs.outputs.values.toSeq)

    markDuplicates(outputFilePaths)

    val moves = outputFilePaths.toList map { case (srcPath, dstPath) =>
      asyncIo.copyAsync(srcPath, dstPath) flatMap { _ =>
        asyncIo.deleteAsync(srcPath)
      }
    }

    Future.sequence(moves)
  }

  /**
    * Happens after everything else runs
    */
  final def afterAll()(implicit ec: ExecutionContext): Future[FinalizationResponse] =
    (workflowDescriptor.finalWorkflowOutputsDir, workflowDescriptor.finalWorkflowOutputsMode) match {
      case (Some(outputsDir), Copy) => copyWorkflowOutputs(outputsDir) map { _ => FinalizationSuccess }
      case (Some(outputsDir), Move) => moveWorkflowOutputs(outputsDir) map { _ => FinalizationSuccess }
      case _ => Future.successful(FinalizationSuccess)
    }
}
