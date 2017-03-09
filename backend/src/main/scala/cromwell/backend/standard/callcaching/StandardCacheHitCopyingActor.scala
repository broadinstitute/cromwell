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
  case object WaitingForCopyResponses extends StandardCacheHitCopyingActorState

  case class StandardCacheHitCopyingActorData(copyCommandsToWaitFor: Set[IoCopyCommand],
                                              copiedJobOutputs: CallOutputs,
                                              copiedDetritus: DetritusMap,
                                              returnCode: Option[Int]
                                             ) {
    def remove(copyCommand: IoCopyCommand) = copy(copyCommandsToWaitFor = copyCommandsToWaitFor filterNot { _ == copyCommand })
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
        case Success((callOutputs, destinationDetritus, allCopyPairs)) =>
          duplicate(allCopyPairs) match {
            case Some(Success(_)) => succeedAndStop(returnCode, callOutputs, destinationDetritus)
            case Some(Failure(failure)) => failAndStop(failure)
            case None =>
              val allCopyCommands = allCopyPairs map { case (source, destination) => copyCommand(source, destination, overwrite = true) }

              allCopyCommands foreach { sendIoCommand(_) }

              goto(WaitingForCopyResponses) using Option(StandardCacheHitCopyingActorData(allCopyCommands, callOutputs, destinationDetritus, returnCode))
          }

        case Failure(failure) => failAndStop(failure)
      }
  }

  when(WaitingForCopyResponses) {
    case Event(IoSuccess(copyCommand: IoCopyCommand, _), Some(data)) =>
      val newData = data.remove(copyCommand)
      if (newData.copyCommandsToWaitFor.isEmpty) succeedAndStop(data.returnCode, data.copiedJobOutputs, data.copiedDetritus)
      else stay() using Option(newData)
    case Event(IoFailure(copyCommand: IoCopyCommand, failure), _) =>
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
    context.parent ! JobSucceededResponse(jobDescriptor.key, returnCode, copiedJobOutputs, Option(detritusMap), Seq.empty)
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

  /**
    * Returns a pair of the list of simpletons with copied paths, and copy commands necessary to perform those copies. 
    */
  private def processSimpletons(wdlValueSimpletons: Seq[WdlValueSimpleton], sourceCallRootPath: Path): Try[(CallOutputs, Set[PathPair])] = Try {
    val (destinationSimpletons, ioCommands): (List[WdlValueSimpleton], Set[PathPair]) = wdlValueSimpletons.toList.foldMap({
      case WdlValueSimpleton(key, wdlFile: WdlFile) =>
        val sourcePath = getPath(wdlFile.value).get
        val destinationPath = PathCopier.getDestinationFilePath(sourceCallRootPath, sourcePath, destinationCallRootPath)

        val destinationSimpleton = WdlValueSimpleton(key, WdlFile(destinationPath.pathAsString))

        List(destinationSimpleton) -> Set(sourcePath -> destinationPath)
      case nonFileSimpleton => (List(nonFileSimpleton), Set.empty[PathPair])
    })

    (WdlValueBuilder.toJobOutputs(jobDescriptor.call.task.outputs, destinationSimpletons), ioCommands)
  }

  /**
    * Returns a pair of the detritus with copied paths, and copy commands necessary to perform those copies. 
    */
  private def processDetritus(sourceJobDetritusFiles: Map[String, String]): Try[(Map[String, Path], Set[PathPair])] = Try {
    val sourceKeys = sourceJobDetritusFiles.keySet
    val destinationKeys = destinationJobDetritusPaths.keySet
    val fileKeys = sourceKeys.intersect(destinationKeys).filterNot(_ == JobPaths.CallRootPathKey)

    val zero = (Map.empty[String, Path], Set.empty[PathPair])

    val (destinationDetritus, ioCommands) = fileKeys.foldLeft(zero)({
      case ((detrituses, commands), detritus) =>
        val sourcePath = getPath(sourceJobDetritusFiles(detritus)).get
        val destinationPath = destinationJobDetritusPaths(detritus)
        
        val newDetrituses = detrituses + (detritus -> destinationPath)
        
        (newDetrituses, commands + ((sourcePath, destinationPath)))
    })
    
    (destinationDetritus + (JobPaths.CallRootPathKey -> destinationCallRootPath), ioCommands)
  }

  override protected def onTimeout(message: Any, to: ActorRef): Unit = {
    val exceptionMessage = message match {
      case copyCommand: IoCopyCommand => s"The Cache hit copying actor timed out waiting for a response to copy ${copyCommand.source.pathAsString} to ${copyCommand.destination.pathAsString}"
      case other => s"The Cache hit copying actor timed out waiting for an unknown I/O operation: $other"
    }

    failAndStop(new TimeoutException(exceptionMessage))
    ()
  }
}
