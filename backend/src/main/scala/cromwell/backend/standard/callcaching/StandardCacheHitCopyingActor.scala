package cromwell.backend.standard.callcaching

import java.util.concurrent.TimeoutException

import akka.actor.{ActorRef, FSM}
import cats.instances.list._
import cats.instances.set._
import cats.instances.tuple._
import cats.syntax.foldable._
import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, JobFailedNonRetryableResponse, JobSucceededResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.backend.standard.callcaching.StandardCacheHitCopyingActor._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}
import cromwell.core._
import cromwell.core.io._
import cromwell.core.logging.JobLogging
import cromwell.core.path.{Path, PathCopier}
import cromwell.core.simpleton.{WdlValueBuilder, WdlValueSimpleton}
import wdl4s.values.WdlFile

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

/**
  * Trait of parameters passed to a StandardCacheHitCopyingActor.
  */
trait StandardCacheHitCopyingActorParams {
  def jobDescriptor: BackendJobDescriptor

  def backendInitializationDataOption: Option[BackendInitializationData]

  def serviceRegistryActor: ActorRef

  def ioActor: ActorRef

  def configurationDescriptor: BackendConfigurationDescriptor
}

/** A default implementation of the cache hit copying params. */
case class DefaultStandardCacheHitCopyingActorParams
(
  override val jobDescriptor: BackendJobDescriptor,
  override val backendInitializationDataOption: Option[BackendInitializationData],
  override val serviceRegistryActor: ActorRef,
  override val ioActor: ActorRef,
  override val configurationDescriptor: BackendConfigurationDescriptor
) extends StandardCacheHitCopyingActorParams

object StandardCacheHitCopyingActor {
  type DetritusMap = Map[String, Path]
  type PathPair = (Path, Path)

  sealed trait StandardCacheHitCopyingActorState
  case object Idle extends StandardCacheHitCopyingActorState
  case object WaitingForIoResponses extends StandardCacheHitCopyingActorState
  case object FailedState extends StandardCacheHitCopyingActorState
  case object WaitingForOnSuccessResponse extends StandardCacheHitCopyingActorState

  case class StandardCacheHitCopyingActorData(commandsToWaitFor: Set[IoCommand[_]],
                                              newJobOutputs: CallOutputs,
                                              newDetritus: DetritusMap,
                                              returnCode: Option[Int]
                                             ) {
    def remove(command: IoCommand[_]) = copy(commandsToWaitFor = commandsToWaitFor - command)
  }
}

class DefaultStandardCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams) extends StandardCacheHitCopyingActor(standardParams) with DefaultIoCommandBuilder

/**
  * Standard implementation of a BackendCacheHitCopyingActor.
  */
abstract class StandardCacheHitCopyingActor(val standardParams: StandardCacheHitCopyingActorParams)
  extends FSM[StandardCacheHitCopyingActorState, Option[StandardCacheHitCopyingActorData]] with JobLogging with StandardCachingActorHelper with IoClientHelper { this: IoCommandBuilder =>

  override lazy val jobDescriptor: BackendJobDescriptor = standardParams.jobDescriptor
  override lazy val backendInitializationDataOption: Option[BackendInitializationData] = standardParams.backendInitializationDataOption
  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor
  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor

  lazy val destinationCallRootPath: Path = jobPaths.callRoot
  lazy val destinationJobDetritusPaths: Map[String, Path] = jobPaths.detritusPaths
  lazy val ioActor = standardParams.ioActor

  startWith(Idle, None)

  context.become(ioReceive orElse receive)

  /** Override this method if you want to provide an alternative way to duplicate files than copying them. */
  protected def duplicate(copyPairs: Set[PathPair]): Option[Try[Unit]] = None

  when(Idle) {
    case Event(CopyOutputsCommand(simpletons, jobDetritus, returnCode), None) =>
      val sourceCallRootPath = lookupSourceCallRootPath(jobDetritus)

      val processed = for {
        (callOutputs, simpletonCopyPairs) <- processSimpletons(simpletons, sourceCallRootPath)
        (destinationDetritus, detritusCopyPairs) <- processDetritus(jobDetritus)
      } yield (callOutputs, destinationDetritus, simpletonCopyPairs ++ detritusCopyPairs)

      processed match {
        case Success((callOutputs, destinationDetritus, allIoCommands)) =>
          duplicate(ioCommandsToCopyPairs(allIoCommands)) match {
            case Some(Success(_)) => succeedAndStop(returnCode, callOutputs, destinationDetritus)
            case Some(Failure(failure)) => failAndStop(failure)
            case None =>

              if (allIoCommands.nonEmpty) {
                allIoCommands foreach { sendIoCommand(_) }
                
                goto(WaitingForIoResponses) using Option(StandardCacheHitCopyingActorData(allIoCommands, callOutputs, destinationDetritus, returnCode))
              } else succeedAndStop(returnCode, callOutputs, destinationDetritus)
          }

        case Failure(failure) => failAndStop(failure)
      }
  }

  when(WaitingForIoResponses) {
    case Event(IoSuccess(command: IoCommand[_], _), Some(data)) =>
      val newData = data.remove(command)
      if (newData.commandsToWaitFor.isEmpty) {
        onSuccessIoCommand(newData) match {
          case Some(successCommand) =>
            sendIoCommand(successCommand)
            goto(WaitingForOnSuccessResponse)
          case None => succeedAndStop(data.returnCode, data.newJobOutputs, data.newDetritus)
        }
      }
      else stay() using Option(newData)
    case Event(IoFailure(_: IoCommand[_], failure), None) =>
      failAndStop(failure)
    case Event(IoFailure(_: IoCommand[_], failure), Some(data)) if data.commandsToWaitFor.nonEmpty =>
      context.parent ! JobFailedNonRetryableResponse(jobDescriptor.key, failure, None)
      // Wait for the other responses to avoid them being sent to dead letter
      goto(FailedState)
    case Event(IoFailure(_: IoCommand[_], failure), Some(_)) =>
      failAndStop(failure)
  }
  
  when(FailedState) {
    case Event(IoFailure(_: IoCommand[_], _), Some(data)) if data.commandsToWaitFor.nonEmpty =>
      stay()
    case Event(IoFailure(_: IoCommand[_], _), Some(_)) =>
      context stop self
      stay()
  }
  
  when(WaitingForOnSuccessResponse) {
    case Event(IoSuccess(_: IoCommand[_], _), Some(data)) =>
      succeedAndStop(data.returnCode, data.newJobOutputs, data.newDetritus)
    case Event(IoFailure(_: IoCommand[_], failure), _) =>
      failAndStop(failure)
  }

  whenUnhandled {
    case Event(AbortJobCommand, _) =>
      abort()
    case Event(unexpected, _) =>
      log.warning(s"Backend cache hit copying actor received an unexpected message: $unexpected in state $stateName")
      stay()
  }

  def succeedAndStop(returnCode: Option[Int], copiedJobOutputs: CallOutputs, detritusMap: DetritusMap) = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), startMetadataKeyValues)
    context.parent ! JobSucceededResponse(jobDescriptor.key, returnCode, copiedJobOutputs, Option(detritusMap), Seq.empty, None)
    context stop self
    stay()
  }

  def failAndStop(failure: Throwable) = {
    context.parent ! JobFailedNonRetryableResponse(jobDescriptor.key, failure, None)
    context stop self
    stay()
  }

  def abort() = {
    log.warning("{}: Abort not supported during cache hit copying", jobTag)
    context.parent ! AbortedResponse(jobDescriptor.key)
    context stop self
    stay()
  }

  private def lookupSourceCallRootPath(sourceJobDetritusFiles: Map[String, String]): Path = {
    sourceJobDetritusFiles.get(JobPaths.CallRootPathKey).map(getPath).get recover {
      case failure =>
        throw new RuntimeException(s"${JobPaths.CallRootPathKey} wasn't found for call ${jobDescriptor.call.fullyQualifiedName}", failure)
    } get
  }
  
  private def ioCommandsToCopyPairs(commands: Set[IoCommand[_]]): Set[PathPair] = commands collect {
    case copyCommand: IoCopyCommand => copyCommand.source -> copyCommand.destination
  }

  /**
    * Returns a pair of the list of simpletons with copied paths, and copy commands necessary to perform those copies. 
    */
  protected def processSimpletons(wdlValueSimpletons: Seq[WdlValueSimpleton], sourceCallRootPath: Path): Try[(CallOutputs, Set[IoCommand[_]])] = Try {
    val (destinationSimpletons, ioCommands): (List[WdlValueSimpleton], Set[IoCommand[_]]) = wdlValueSimpletons.toList.foldMap({
      case WdlValueSimpleton(key, wdlFile: WdlFile) =>
        val sourcePath = getPath(wdlFile.value).get
        val destinationPath = PathCopier.getDestinationFilePath(sourceCallRootPath, sourcePath, destinationCallRootPath)

        val destinationSimpleton = WdlValueSimpleton(key, WdlFile(destinationPath.pathAsString))

        List(destinationSimpleton) -> Set(copyCommand(sourcePath, destinationPath, overwrite = true))
      case nonFileSimpleton => (List(nonFileSimpleton), Set.empty[IoCommand[_]])
    })

    (WdlValueBuilder.toJobOutputs(jobDescriptor.call.task.outputs, destinationSimpletons), ioCommands)
  }

  protected final def detritusFileKeys(sourceJobDetritusFiles: Map[String, String]) = {
    val sourceKeys = sourceJobDetritusFiles.keySet
    val destinationKeys = destinationJobDetritusPaths.keySet
    sourceKeys.intersect(destinationKeys).filterNot(_ == JobPaths.CallRootPathKey)
  }
  
  /**
    * Returns a pair of the detritus with copied paths, and copy commands necessary to perform those copies. 
    */
  protected def processDetritus(sourceJobDetritusFiles: Map[String, String]): Try[(Map[String, Path], Set[IoCommand[_]])] = Try {
    val fileKeys = detritusFileKeys(sourceJobDetritusFiles)

    val zero = (Map.empty[String, Path], Set.empty[IoCommand[_]])

    val (destinationDetritus, ioCommands) = fileKeys.foldLeft(zero)({
      case ((detrituses, commands), detritus) =>
        val sourcePath = getPath(sourceJobDetritusFiles(detritus)).get
        val destinationPath = destinationJobDetritusPaths(detritus)
        
        val newDetrituses = detrituses + (detritus -> destinationPath)
        
        (newDetrituses, commands + copyCommand(sourcePath, destinationPath, overwrite = true))
    })
    
    (destinationDetritus + (JobPaths.CallRootPathKey -> destinationCallRootPath), ioCommands)
  }
  
  protected def onSuccessIoCommand(data: StandardCacheHitCopyingActorData): Option[IoCommand[_]] = None

  override protected def onTimeout(message: Any, to: ActorRef): Unit = {
    val exceptionMessage = message match {
      case copyCommand: IoCopyCommand => s"The Cache hit copying actor timed out waiting for a response to copy ${copyCommand.source.pathAsString} to ${copyCommand.destination.pathAsString}"
      case other => s"The Cache hit copying actor timed out waiting for an unknown I/O operation: $other"
    }

    failAndStop(new TimeoutException(exceptionMessage))
    ()
  }
}
