package cromwell.engine.workflow.lifecycle.deletion

import java.io.FileNotFoundException

import akka.actor.{ActorRef, LoggingFSM, Props}
import cromwell.core.io._
import cromwell.core.path.{Path, PathBuilder, PathFactory}
import cromwell.core.{RootWorkflowId, WorkflowId, WorkflowMetadataKeys}
import cromwell.engine.io.IoAttempts.EnhancedCromwellIoException
import cromwell.engine.workflow.lifecycle.deletion.DeleteWorkflowFilesActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching._
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import cromwell.services.EngineServicesStore
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.impl.FileDeletionStatus
import cromwell.services.metadata.impl.FileDeletionStatus.{Failed, InProgress, Succeeded}
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.apache.commons.lang3.exception.ExceptionUtils
import wom.values.{WomSingleFile, WomValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

//noinspection DuplicatedCode
class DeleteWorkflowFilesActor(rootWorkflowId: RootWorkflowId,
                               rootAndSubworkflowIds: Set[WorkflowId],
                               workflowFinalOutputs: Set[WomValue],
                               workflowAllOutputs: Set[WomValue],
                               pathBuilders: List[PathBuilder],
                               serviceRegistryActor: ActorRef,
                               ioActorRef: ActorRef) extends LoggingFSM[DeleteWorkflowFilesActorState, DeleteWorkflowFilesActorStateData] with IoClientHelper {

  implicit val ec: ExecutionContext = context.dispatcher

  val gcsCommandBuilder = GcsBatchCommandBuilder
  val asyncIO = new AsyncIo(ioActorRef, gcsCommandBuilder)
  val callCache = new CallCache(EngineServicesStore.engineDatabaseInterface)

  startWith(Pending, NoData)

  when(Pending) {
    case Event(StartWorkflowFilesDeletion, NoData) =>
      val intermediateOutputs = gatherIntermediateOutputFiles(workflowAllOutputs, workflowFinalOutputs)
      if (intermediateOutputs.nonEmpty) {
        self ! DeleteFiles
        goto(DeleteIntermediateFiles) using DeletingIntermediateFilesData(intermediateOutputs)
      }
      else {
        log.info(s"Root workflow ${rootWorkflowId.id} does not have any intermediate output files to delete.")
        respondAndStop(Nil, Nil, Nil)
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
    case Event(InvalidateCallCache, _) =>
      fetchCallCacheEntries(callCache) onComplete {
        case Failure(throwable) => self ! FailedRetrieveCallCacheIds(throwable)
        case Success(ids) => self ! RetrievedCallCacheIds(ids)
      }
      stay()
    case Event(FailedRetrieveCallCacheIds(throwable), data: InvalidateCallCacheData) =>
      log.error(
        throwable,
        s"Unable to retrieve call cache values for for root and subworkflows: ${rootAndSubworkflowIds.mkString(",")}."
      )
      respondAndStop(data.deleteErrors, data.filesNotFound, List.empty)
    case Event(RetrievedCallCacheIds(callCacheEntryIds), data: InvalidateCallCacheData) if callCacheEntryIds.isEmpty =>
      log.error(
        s"No call cache entries found to invalidate for root and subworkflows: ${rootAndSubworkflowIds.mkString(",")}."
      )
      respondAndStop(data.deleteErrors, data.filesNotFound, List.empty)
    case Event(RetrievedCallCacheIds(callCacheEntryIds), data: InvalidateCallCacheData) =>
      callCacheEntryIds.map(id => context.actorOf(CallCacheInvalidateActor.props(callCache, CallCachingEntryId(id))))
      goto(WaitingForInvalidateCCResponses) using
        WaitingForInvalidateCCResponsesData(callCacheEntryIds, data.deleteErrors, data.filesNotFound)
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

  private def toPath(womSingleFile: WomSingleFile): Option[Path] = {
    Try(PathFactory.buildPath(womSingleFile.valueString, pathBuilders)).toOption
  }

  private def getWomSingleFiles(womValue: WomValue): Seq[WomSingleFile] = {
    womValue.collectAsSeq({ case womSingleFile: WomSingleFile => womSingleFile })
  }

  def gatherIntermediateOutputFiles(allOutputs: Set[WomValue], finalOutputs: Set[WomValue]): Set[Path] = {
    val allOutputFiles = allOutputs.flatMap(getWomSingleFiles)
    val finalOutputFiles = finalOutputs.flatMap(getWomSingleFiles)
    allOutputFiles.diff(finalOutputFiles).flatMap(toPath)
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

  //@formatter:off
  // Commands
  sealed trait DeleteWorkflowFilesActorMessage
  object StartWorkflowFilesDeletion extends DeleteWorkflowFilesActorMessage
  object DeleteFiles extends DeleteWorkflowFilesActorMessage
  object InvalidateCallCache extends DeleteWorkflowFilesActorMessage
  case class RetrievedCallCacheIds(ids: Set[Int]) extends DeleteWorkflowFilesActorMessage
  case class FailedRetrieveCallCacheIds(throwable: Throwable) extends DeleteWorkflowFilesActorMessage

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
  //@formatter:on

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
            workflowFinalOutputs: Set[WomValue],
            workflowAllOutputs: Set[WomValue],
            pathBuilders: List[PathBuilder],
            serviceRegistryActor: ActorRef,
            ioActor: ActorRef): Props = {
    Props(new DeleteWorkflowFilesActor(rootWorkflowId, rootAndSubworkflowIds, workflowFinalOutputs, workflowAllOutputs, pathBuilders, serviceRegistryActor, ioActor))
  }
}
