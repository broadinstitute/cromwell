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

  // TODO: this mechanism here is very close to the one in CallCacheHashingJobActorData
  // Abstracting it might be valuable
  /**
    * The head subset of commandsToWaitFor is sent to the IoActor as a bulk.
    * When a response comes back, the corresponding command is removed from the head set.
    * When the head set is empty, it is removed and the next subset is sent, until there is no subset left.
    * If at any point a response comes back as a failure. Other responses for the current set will be awaited for 
    * but subsequent sets will not be sent and the actor will send back a failure message.
    */
  case class StandardCacheHitCopyingActorData(commandsToWaitFor: List[Set[IoCommand[_]]],
                                              newJobOutputs: CallOutputs,
                                              newDetritus: DetritusMap,
                                              returnCode: Option[Int]
                                             ) {

    /**
      * Removes the command from commandsToWaitFor
      * returns a pair of the new state data and CommandSetState giving information about what to do next
      */
    def commandComplete(command: IoCommand[_]): (StandardCacheHitCopyingActorData, CommandSetState) = commandsToWaitFor match {
      // If everything was already done send back current data and AllCommandsDone
      case Nil => (this, AllCommandsDone)
      case lastSubset :: Nil =>
        val updatedSubset = lastSubset - command
        // If the last subset is now empty, we're done
        if (updatedSubset.isEmpty) (this.copy(commandsToWaitFor = List.empty), AllCommandsDone)
        // otherwise update commandsToWaitFor and keep waiting
        else (this.copy(commandsToWaitFor = List(updatedSubset)), StillWaiting)
      case currentSubset :: otherSubsets =>
        val updatedSubset = currentSubset - command
        // This subset is done but there are other ones, remove it from commandsToWaitFor and return the next round of commands
        if (updatedSubset.isEmpty) (this.copy(commandsToWaitFor = otherSubsets), NextSubSet(otherSubsets.head))
        // otherwise update the head susbset and keep waiting  
        else (this.copy(commandsToWaitFor = List(updatedSubset) ++ otherSubsets), StillWaiting)
    }
  }

  // Internal ADT to keep track of command set states
  private[callcaching] sealed trait CommandSetState
  private[callcaching] case object StillWaiting extends CommandSetState
  private[callcaching] case object AllCommandsDone extends CommandSetState
  private[callcaching] case class NextSubSet(commands: Set[IoCommand[_]]) extends CommandSetState
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

      // Try to make a Path of the callRootPath from the detritus
      lookupSourceCallRootPath(jobDetritus) match {
        case Success(sourceCallRootPath) =>
          
          // process simpletons and detritus to get updated paths and corresponding IoCommands
          val processed = for {
            (destinationCallOutputs, simpletonIoCommands) <- processSimpletons(simpletons, sourceCallRootPath)
            (destinationDetritus, detritusIoCommands) <- processDetritus(jobDetritus)
          } yield (destinationCallOutputs, destinationDetritus, simpletonIoCommands ++ detritusIoCommands)

          processed match {
            case Success((destinationCallOutputs, destinationDetritus, detritusAndOutputsIoCommands)) =>
              duplicate(ioCommandsToCopyPairs(detritusAndOutputsIoCommands)) match {
                  // Use the duplicate override if exists
                case Some(Success(_)) => succeedAndStop(returnCode, destinationCallOutputs, destinationDetritus)
                case Some(Failure(failure)) => failAndStop(failure)
                  // Otherwise send the first round of IoCommands (file outputs and detritus) if any
                case None if detritusAndOutputsIoCommands.nonEmpty =>
                    detritusAndOutputsIoCommands foreach { sendIoCommand(_) }

                    // Add potential additional commands to the list
                  val additionalCommands = additionalIoCommands(sourceCallRootPath, simpletons, destinationCallOutputs, jobDetritus, destinationDetritus)
                  val allCommands = List(detritusAndOutputsIoCommands) ++ additionalCommands

                    goto(WaitingForIoResponses) using Option(StandardCacheHitCopyingActorData(allCommands, destinationCallOutputs, destinationDetritus, returnCode))
                case _ => succeedAndStop(returnCode, destinationCallOutputs, destinationDetritus)
              }

            case Failure(failure) => failAndStop(failure)
          }

        case Failure(failure) => failAndStop(failure)
      }
  }

  when(WaitingForIoResponses) {
    case Event(IoSuccess(command: IoCommand[_], _), Some(data)) =>
      val (newData, commandState) = data.commandComplete(command)

      commandState match {
        case StillWaiting => stay() using Option(newData)
        case AllCommandsDone => succeedAndStop(newData.returnCode, newData.newJobOutputs, newData.newDetritus)
        case NextSubSet(commands) =>
          commands foreach { sendIoCommand(_) }
          stay() using Option(newData)
      }
    case Event(IoFailure(command: IoCommand[_], failure), Some(data)) =>
      // any failure is fatal
      context.parent ! JobFailedNonRetryableResponse(jobDescriptor.key, failure, None)

      val (newData, commandState) = data.commandComplete(command)

      commandState match {
        // If we're still waiting for some responses, go to failed state
        case StillWaiting => goto(FailedState) using Option(newData)
        // Otherwise we're done
        case _ =>
          context stop self
          stay()
      }
    // Should not be possible
    case Event(IoFailure(_: IoCommand[_], failure), None) => failAndStop(failure)
  }

  when(FailedState) {
    // At this point success or failure doesn't matter, we've already failed this hit
    case Event(response: IoAck[_], Some(data)) =>
      val (newData, commandState) = data.commandComplete(response.command)
      commandState match {
        // If we're still waiting for some responses, stay
        case StillWaiting => stay() using Option(newData)
        // Otherwise we're done
        case _ =>
          context stop self
          stay()
      }
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

  protected def lookupSourceCallRootPath(sourceJobDetritusFiles: Map[String, String]): Try[Path] = {
    sourceJobDetritusFiles.get(JobPaths.CallRootPathKey) match {
      case Some(source) => getPath(source)
      case None => Failure(new RuntimeException(s"${JobPaths.CallRootPathKey} wasn't found for call ${jobDescriptor.call.fullyQualifiedName}"))
    }
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

  /**
    * Returns the file (and ONLY the file detritus) intersection between the cache hit and this call.
    */
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

  /**
    * Additional IoCommands that will be sent after (and only after) output and detritus commands complete successfully.
    * See StandardCacheHitCopyingActorData
    */
  protected def additionalIoCommands(sourceCallRootPath: Path,
                                     originalSimpletons: Seq[WdlValueSimpleton],
                                     newOutputs: CallOutputs,
                                     originalDetritus:  Map[String, String],
                                     newDetritus: Map[String, Path]): List[Set[IoCommand[_]]] = List.empty

  override protected def onTimeout(message: Any, to: ActorRef): Unit = {
    val exceptionMessage = message match {
      case copyCommand: IoCopyCommand => s"The Cache hit copying actor timed out waiting for a response to copy ${copyCommand.source.pathAsString} to ${copyCommand.destination.pathAsString}"
      case touchCommand: IoTouchCommand => s"The Cache hit copying actor timed out waiting for a response to touch ${touchCommand.file.pathAsString}"
      case other => s"The Cache hit copying actor timed out waiting for an unknown I/O operation: $other"
    }

    failAndStop(new TimeoutException(exceptionMessage))
    ()
  }
}
