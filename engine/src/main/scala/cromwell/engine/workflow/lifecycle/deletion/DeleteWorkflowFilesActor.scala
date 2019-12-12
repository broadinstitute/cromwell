package cromwell.engine.workflow.lifecycle.deletion

import java.io.FileNotFoundException

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.io._
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder, PathFactory}
import cromwell.core.{CallOutputs, RootWorkflowId, WorkflowId, WorkflowMetadataKeys}
import cromwell.engine.io.IoAttempts.EnhancedCromwellIoException
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCache, CallCacheInvalidateActor, CallCacheInvalidatedFailure, CallCacheInvalidatedSuccess, CallCachingEntryId}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import cromwell.services.EngineServicesStore
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.impl.FileDeletionStatus
import cromwell.services.metadata.impl.FileDeletionStatus.{Failed, InProgress, Succeeded}
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.lang3.exception.ExceptionUtils
import wom.graph.GraphNodePort.OutputPort
import wom.types.WomSingleFileType
import wom.values.{WomArray, WomSingleFile, WomValue}

import scala.concurrent.{Await, ExecutionContext, Future}
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

class DeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId,
                               rootAndSubworkflowIds: Set[WorkflowId],
                               workflowFinalOutputs: CallOutputs,
                               workflowAllOutputs: CallOutputs,
                               pathBuilders: List[PathBuilder],
                               serviceRegistryActor: ActorRef,
                               ioActorRef: ActorRef) extends LoggingFSM[DeleteWorkflowFilesActorState, DeleteWorkflowFilesActorStateData] with IoClientHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  /*
  building a path for a string output such as 'Hello World' satisfies the conditions of a
  being a relative file path for a DefaultPath which we don't want, so remove DefaultPathBuilder
  (this is because we only want cloud backend paths like gs://, s3://, etc.)
   */
  val pathBuildersWithoutDefault = pathBuilders.filterNot(_ == DefaultPathBuilder)
  val gcsCommandBuilder = GcsBatchCommandBuilder
  val asyncIO = new AsyncIo(ioActorRef, gcsCommandBuilder)

  startWith(Pending, NoData)

  when(Pending) {
    case Event(StartWorkflowFilesDeletion, NoData) =>
      val intermediateOutputs = gatherIntermediateOutputFiles(workflowAllOutputs.outputs, workflowFinalOutputs.outputs)
      if (intermediateOutputs.nonEmpty) {
        self ! DeleteFiles
        goto(DeleteIntermediateFiles) using DeletingIntermediateFilesData(intermediateOutputs)
      }
      else {
        log.info(s"Root workflow ${rootWorkflowId.id} does not have any intermediate output files to delete.")
        stopSelf()
      }
  }

  when(DeleteIntermediateFiles) {
    case Event(DeleteFiles, DeletingIntermediateFilesData(intermediateFiles)) =>
      // update deletion status in metadata
      val deletionInProgressEvent = metadataEventForDeletionStatus(InProgress)
      serviceRegistryActor ! PutMetadataAction(deletionInProgressEvent)

      // send delete IoCommand for each file to ioActor
      val deleteCommands: Set[IoDeleteCommand] = intermediateFiles.map(gcsCommandBuilder.deleteCommand(_, swallowIoExceptions = false))
      deleteCommands foreach sendIoCommand

      goto(WaitingForIoResponses) using WaitingForIoResponsesData(deleteCommands)
  }

  when(WaitingForIoResponses) {
    case Event(IoSuccess(command: IoDeleteCommand, _), data: WaitingForIoResponsesData) =>
      val (newData: WaitingForIoResponsesData, commandState) = data.commandComplete(command)
      commandState match {
        case StillWaiting => stay() using newData
        case AllCommandsDone =>
          // once deletion is complete, invalidate call cache entries
          self ! InvalidateCallCache
          goto(InvalidatingCallCache) using InvalidateCallCacheData(newData.deleteErrors, newData.filesNotFound)
      }
    case Event(IoFailure(command: IoDeleteCommand, error: Throwable), data: WaitingForIoResponsesData) =>
      val (newData: WaitingForIoResponsesData, commandState) = data.commandComplete(command)

      /*
        FileNotFound exception might be because of:
          - file does not exist because of some reason (maybe user deleted it?)
          - it was deleted by Cromwell, but before the status of the workflow could be updated by WorkflowActor,
            Cromwell was restarted. So upon restart the deletion process is ran again entirely and hence now previously
            deleted files will not exist.
         In both these cases, we consider the deletion process a success, but warn the users of such files not found.
       */
      val newDataWithErrorUpdates = error match {
        case EnhancedCromwellIoException(_, _: FileNotFoundException) => newData.copy(filesNotFound = newData.filesNotFound :+ command.file)
        case _ => newData.copy(deleteErrors = newData.deleteErrors :+ error)
      }
      commandState match {
        case StillWaiting => stay() using newDataWithErrorUpdates
        case AllCommandsDone =>
          // once deletion is complete, invalidate call cache entries
          self ! InvalidateCallCache
          goto(InvalidatingCallCache) using InvalidateCallCacheData(newDataWithErrorUpdates.deleteErrors, newDataWithErrorUpdates.filesNotFound)
      }
  }

  when(InvalidatingCallCache) {
    case Event(InvalidateCallCache, data: InvalidateCallCacheData) => {
      val callCache = new CallCache(EngineServicesStore.engineDatabaseInterface)
      val callCacheEntryIds = Await.result(fetchCallCacheEntries(callCache), 20.seconds)

      if (callCacheEntryIds.isEmpty) {
        log.error(s"No call cache entries found to invalidate for root and subworkflows: ${rootAndSubworkflowIds.mkString(",")}.")
        respondAndStop(data.deleteErrors, data.filesNotFound, List.empty)
      } else {
        // create CallCacheInvalidateActor for each call cache entry id and wait for its response
        callCacheEntryIds.map(id => context.actorOf(CallCacheInvalidateActor.props(callCache, CallCachingEntryId(id))))
        goto(WaitingForInvalidateCCResponses) using WaitingForInvalidateCCResponsesData(callCacheEntryIds, data.deleteErrors, data.filesNotFound)
      }
    }
  }

  when(WaitingForInvalidateCCResponses) {
    case Event(CallCacheInvalidatedSuccess(cacheId, _), data: WaitingForInvalidateCCResponsesData) =>
      val (newData: WaitingForInvalidateCCResponsesData, invalidateState) = data.commandComplete(cacheId.id)
      invalidateState match {
        case StillWaiting => stay() using newData
        case AllCommandsDone => respondAndStop(newData.deleteErrors, newData.filesNotFound, newData.callCacheInvalidationErrors)
      }
    case Event(CallCacheInvalidatedFailure(cacheId, error), data: WaitingForInvalidateCCResponsesData) =>
      val (newData: WaitingForInvalidateCCResponsesData, invalidateState) = data.commandComplete(cacheId.id)
      val updatedDataWithError = newData.copy(callCacheInvalidationErrors = newData.callCacheInvalidationErrors :+ error)
      invalidateState match {
        case StillWaiting => stay() using updatedDataWithError
        case AllCommandsDone => respondAndStop(newData.deleteErrors, newData.filesNotFound, updatedDataWithError.callCacheInvalidationErrors)
      }
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


  private def respondAndStop(errors: List[Throwable], filesNotFound: List[Path], callCacheInvalidationErrors: List[Throwable]) = {
    val (metadataEvent, response) =
      if (errors.isEmpty) (metadataEventForDeletionStatus(Succeeded), DeleteWorkflowFilesSucceededResponse(filesNotFound, callCacheInvalidationErrors))
      else (metadataEventForDeletionStatus(Failed), DeleteWorkflowFilesFailedResponse(errors, filesNotFound, callCacheInvalidationErrors))

    serviceRegistryActor ! PutMetadataAction(metadataEvent)
    context.parent ! response
    stopSelf()
  }


  private def metadataEventForDeletionStatus(status: FileDeletionStatus): MetadataEvent = {
    val key = MetadataKey(rootWorkflowId, None, WorkflowMetadataKeys.FileDeletionStatus)
    val value = MetadataValue(FileDeletionStatus.toDatabaseValue(status))

    MetadataEvent(key, value)
  }


  private def fetchCallCacheEntries(callCache: CallCache): Future[Set[Int]] = {
    val callCacheEntryIdsFuture = rootAndSubworkflowIds.map(x => callCache.callCacheEntryIdsForWorkflowId(x.toString)).map { f =>
      f.map { Success(_) }.recover { case t => Failure(t) }}

    Future.sequence(callCacheEntryIdsFuture).map { _.flatMap {
      case Success(callCacheEntryIds) =>
        Option(callCacheEntryIds)
      case Failure(e) =>
        log.error(s"Failed to fetch call cache entry ids for workflow. Error: ${ExceptionUtils.getMessage(e)}")
        None
    }.flatten}
  }


  def gatherIntermediateOutputFiles(allOutputs: Map[OutputPort, WomValue], finalOutputs: Map[OutputPort, WomValue]): Set[Path] = {
    val isFinalOutputsNonEmpty = finalOutputs.nonEmpty

    def existsInFinalOutputs(value: String): Option[Boolean] = {
      finalOutputs.map{ case (_, output) => value.equals(output.valueString)}.reduceOption(_ || _)
    }

    def getFilePathForSingleFile(value: String): Option[Path] = {
      //if workflow has final outputs check if this output value does not exist in final outputs
      (isFinalOutputsNonEmpty, existsInFinalOutputs(value)) match {
        case (true, Some(true)) => None
        case (true, Some(false)) | (true, None) | (false, _) =>
          // eliminate outputs which are not files
          Try(PathFactory.buildPath(value, pathBuildersWithoutDefault)).toOption
      }
    }

    def getFilePath(value: WomValue): Seq[Path] = {
      val seqOfPaths = value match {
        case v: WomArray => v.value.map(a => getFilePathForSingleFile(a.valueString))
        case s => Seq(getFilePathForSingleFile(s.valueString))
      }
      seqOfPaths.flatten
    }

    def checkOutputIsFile(value: WomValue): Boolean = {
      value match {
        case _: WomSingleFile => true
        case w: WomArray => w.arrayType.memberType.isCoerceableFrom(WomSingleFileType)
        case _ => false
      }
    }

    allOutputs.collect {
      case (_, output) if checkOutputIsFile(output) => getFilePath(output)
    }.flatten.toSet
  }


  override def ioActor: ActorRef = ioActorRef


  override protected def onTimeout(message: Any, to: ActorRef): Unit = {
    message match {
      case delete: IoDeleteCommand => log.error(s"The DeleteWorkflowFilesActor for root workflow $rootWorkflowId timed out " +
        s"waiting for a response for deleting file ${delete.file}.")
      case other => log.error(s"The DeleteWorkflowFilesActor for root workflow $rootWorkflowId timed out " +
        s"waiting for a response for unknown operation: $other.")
    }
  }
}


object DeleteWorkflowFilesActor {

  // Commands
  sealed trait DeleteWorkflowFilesActorMessage
  object StartWorkflowFilesDeletion extends DeleteWorkflowFilesActorMessage
  object DeleteFiles extends DeleteWorkflowFilesActorMessage
  object InvalidateCallCache extends DeleteWorkflowFilesActorMessage

  // Actor states
  sealed trait DeleteWorkflowFilesActorState
  object Pending extends DeleteWorkflowFilesActorState
  object DeleteIntermediateFiles extends DeleteWorkflowFilesActorState
  object WaitingForIoResponses extends DeleteWorkflowFilesActorState
  object InvalidatingCallCache extends DeleteWorkflowFilesActorState
  object WaitingForInvalidateCCResponses extends DeleteWorkflowFilesActorState

  // State data
  sealed trait DeleteWorkflowFilesActorStateData
  object NoData extends DeleteWorkflowFilesActorStateData
  case class DeletingIntermediateFilesData(intermediateFiles: Set[Path]) extends DeleteWorkflowFilesActorStateData
  case class InvalidateCallCacheData(deleteErrors: List[Throwable], filesNotFound: List[Path]) extends DeleteWorkflowFilesActorStateData

  abstract class WaitingForResponseFromActorData[A](commandsToWaitFor: Set[A]) {

    def assertionFailureMsg(expectedSize: Int, requiredSize: Int): String

    def setCommandsToWaitFor(updatedCommandsToWaitFor: Set[A]): WaitingForResponseFromActorData[A]

    def commandComplete(command: A): (WaitingForResponseFromActorData[A], WaitingForResponseState) = {
      if (commandsToWaitFor.isEmpty) (this, AllCommandsDone)
      else {
        val updatedCommandsSet = commandsToWaitFor - command

        val expectedCommandSetSize = updatedCommandsSet.size
        val requiredCommandSetSize = commandsToWaitFor.size - 1
        require(expectedCommandSetSize == requiredCommandSetSize, assertionFailureMsg(expectedCommandSetSize, requiredCommandSetSize))

        if (updatedCommandsSet.isEmpty) (setCommandsToWaitFor(Set.empty), AllCommandsDone)
        else (setCommandsToWaitFor(updatedCommandsSet), StillWaiting)
      }
    }
  }

  case class WaitingForIoResponsesData(commandsToWaitFor: Set[IoDeleteCommand],
                                       deleteErrors: List[Throwable] = List.empty,
                                       filesNotFound: List[Path] = List.empty)
    extends WaitingForResponseFromActorData[IoDeleteCommand](commandsToWaitFor) with DeleteWorkflowFilesActorStateData {

    override def assertionFailureMsg(expectedSize: Int, requiredSize: Int): String = {
      s"Found updated command set size as $expectedSize instead of $requiredSize. The updated set of commands that " +
        s"DeleteWorkflowFilesActor has to wait for should be 1 less after removing a completed command."
    }

    override def setCommandsToWaitFor(updatedCommandsToWaitFor: Set[IoDeleteCommand]): WaitingForResponseFromActorData[IoDeleteCommand] = {
      this.copy(commandsToWaitFor = updatedCommandsToWaitFor)
    }
  }

  case class WaitingForInvalidateCCResponsesData(commandsToWaitFor: Set[Int],
                                                 deleteErrors: List[Throwable],
                                                 filesNotFound: List[Path],
                                                 callCacheInvalidationErrors: List[Throwable] = List.empty)
    extends WaitingForResponseFromActorData[Int](commandsToWaitFor) with DeleteWorkflowFilesActorStateData {

    override def assertionFailureMsg(expectedSize: Int, requiredSize: Int): String = {
      s"Found updated call cache entries set size as $expectedSize instead of $requiredSize. The updated set of call cache entries" +
        s" that DeleteWorkflowFilesActor has to wait for should be 1 less after a call cache entry is invalidated."
    }

    override def setCommandsToWaitFor(updatedCommandsToWaitFor: Set[Int]): WaitingForResponseFromActorData[Int] = {
      this.copy(commandsToWaitFor = updatedCommandsToWaitFor)
    }
  }

  // Responses
  sealed trait DeleteWorkflowFilesResponse
  case class DeleteWorkflowFilesSucceededResponse(filesNotFound: List[Path], callCacheInvalidationErrors: List[Throwable]) extends DeleteWorkflowFilesResponse
  case class DeleteWorkflowFilesFailedResponse(errors: List[Throwable], filesNotFound: List[Path], callCacheInvalidationErrors: List[Throwable]) extends DeleteWorkflowFilesResponse

  // internal state to keep track of deletion of files and call cache invalidation
  sealed trait WaitingForResponseState
  private[deletion] case object StillWaiting extends WaitingForResponseState
  private[deletion] case object AllCommandsDone extends WaitingForResponseState


  def props(rootWorkflowId: RootWorkflowId,
            rootAndSubworkflowIds: Set[WorkflowId],
            workflowFinalOutputs: CallOutputs,
            workflowAllOutputs: CallOutputs,
            pathBuilders: List[PathBuilder],
            serviceRegistryActor: ActorRef,
            ioActor: ActorRef): Props = {
    Props(new DeleteWorkflowFilesActor(rootWorkflowId, rootAndSubworkflowIds, workflowFinalOutputs, workflowAllOutputs, pathBuilders, serviceRegistryActor, ioActor))
  }
}
