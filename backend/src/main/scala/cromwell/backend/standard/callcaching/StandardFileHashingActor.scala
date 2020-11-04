package cromwell.backend.standard.callcaching

import java.util.concurrent.TimeoutException

import akka.actor.{Actor, ActorLogging, ActorRef, Timers}
import akka.event.LoggingAdapter
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.backend.standard.callcaching.RootWorkflowFileHashCacheActor.IoHashCommandWithContext
import cromwell.backend.standard.callcaching.StandardFileHashingActor._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}
import cromwell.core.JobKey
import cromwell.core.callcaching._
import cromwell.core.io._
import cromwell.core.logging.JobLogging
import wom.values.WomFile

import scala.util.{Failure, Success, Try}

/**
  * Trait of parameters passed to a StandardCacheHitCopyingActor.
  */
trait StandardFileHashingActorParams {
  def jobDescriptor: BackendJobDescriptor

  def backendInitializationDataOption: Option[BackendInitializationData]

  def serviceRegistryActor: ActorRef

  def ioActor: ActorRef

  def configurationDescriptor: BackendConfigurationDescriptor

  def fileHashCachingActor: Option[ActorRef]
}

/** A default implementation of the cache hit copying params. */
case class DefaultStandardFileHashingActorParams
(
  override val jobDescriptor: BackendJobDescriptor,
  override val backendInitializationDataOption: Option[BackendInitializationData],
  override val serviceRegistryActor: ActorRef,
  override val ioActor: ActorRef,
  override val configurationDescriptor: BackendConfigurationDescriptor,
  override val fileHashCachingActor: Option[ActorRef]
) extends StandardFileHashingActorParams

case class FileHashContext(hashKey: HashKey, file: String)

class DefaultStandardFileHashingActor(standardParams: StandardFileHashingActorParams) extends StandardFileHashingActor(standardParams) {
  override val ioCommandBuilder: IoCommandBuilder = DefaultIoCommandBuilder
}

object StandardFileHashingActor {
  case class FileHashingFunction(work: (SingleFileHashRequest, LoggingAdapter) => Try[String])

  sealed trait BackendSpecificHasherCommand { def jobKey: JobKey }
  final case class SingleFileHashRequest(jobKey: JobKey, hashKey: HashKey, file: WomFile, initializationData: Option[BackendInitializationData]) extends BackendSpecificHasherCommand

  sealed trait BackendSpecificHasherResponse extends SuccessfulHashResultMessage
  case class FileHashResponse(hashResult: HashResult) extends BackendSpecificHasherResponse { override def hashes = Set(hashResult) }
}

abstract class StandardFileHashingActor(standardParams: StandardFileHashingActorParams)
  extends Actor
    with ActorLogging
    with JobLogging
    with IoClientHelper
    with StandardCachingActorHelper
    with Timers {
  override lazy val ioActor: ActorRef = standardParams.ioActor
  override lazy val jobDescriptor: BackendJobDescriptor = standardParams.jobDescriptor
  override lazy val backendInitializationDataOption: Option[BackendInitializationData] = standardParams.backendInitializationDataOption
  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor
  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor

  protected def ioCommandBuilder: IoCommandBuilder = DefaultIoCommandBuilder

  def customHashStrategy(fileRequest: SingleFileHashRequest): Option[Try[String]] = None
  
  def fileHashingReceive: Receive = {
    // Hash Request
    case fileRequest: SingleFileHashRequest =>
      customHashStrategy(fileRequest) match {
        case Some(Success(result)) => context.parent ! FileHashResponse(HashResult(fileRequest.hashKey, HashValue(result)))
        case Some(Failure(failure)) => context.parent ! HashingFailedMessage(fileRequest.file.value, failure)
        case None => asyncHashing(fileRequest, context.parent)
      }

    // Hash Success
    case (fileHashRequest: FileHashContext, IoSuccess(_, result: String)) =>
      context.parent ! FileHashResponse(HashResult(fileHashRequest.hashKey, HashValue(result)))

    case (fileHashRequest: FileHashContext, IoSuccess(_, other)) =>
      context.parent ! HashingFailedMessage(
        fileHashRequest.file,
        new Exception(s"Hash function supposedly succeeded but responded with '$other' instead of a string hash"),
      )

    // Hash Failure
    case (fileHashRequest: FileHashContext, IoFailAck(_, failure: Throwable)) =>
      context.parent ! HashingFailedMessage(fileHashRequest.file, failure)

    case other =>
      log.warning(s"Async File hashing actor received unexpected message: $other")
  }

  def asyncHashing(fileRequest: SingleFileHashRequest, replyTo: ActorRef): Unit = {
    val fileAsString = fileRequest.file.value
    val ioHashCommandTry = for {
      gcsPath <- getPath(fileAsString)
      command <- ioCommandBuilder.hashCommand(gcsPath)
    } yield command
    lazy val fileHashContext = FileHashContext(fileRequest.hashKey, fileRequest.file.value)

    ioHashCommandTry match {
      case Success(ioHashCommand) =>
        // If there is a file hash caching actor then forward the IoHashCommandWithContext to that. If there isn't a
        // file hash caching actor then this actor should send the IO command itself.
        standardParams.fileHashCachingActor match {
          case Some(cacheActor) => cacheActor ! IoHashCommandWithContext(ioHashCommand, fileHashContext)
          case None => sendIoCommandWithContext(ioHashCommand, fileHashContext)
        }
      case Failure(failure) => replyTo ! HashingFailedMessage(fileAsString, failure)
    }
  }

  override def receive: Receive = ioReceive orElse fileHashingReceive 

  override protected def onTimeout(message: Any, to: ActorRef): Unit = {
    message match {
      case (_, ioHashCommand: IoHashCommand) =>
        val fileAsString = ioHashCommand.file.pathAsString
        context.parent !
          HashingFailedMessage(fileAsString, new TimeoutException(s"Hashing request timed out for: $fileAsString"))
      case other =>
        // This should never happen... but at least send _something_ before this actor goes silent.
        log.warning(s"Async File hashing actor received unexpected timeout message: $other")
        context.parent ! HashingServiceUnvailable
    }
  }
}
