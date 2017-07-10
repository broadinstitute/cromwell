package cromwell.backend.standard.callcaching

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.event.LoggingAdapter
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.backend.standard.callcaching.StandardFileHashingActor.{FileHashResponse, SingleFileHashRequest}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}
import cromwell.core.JobKey
import cromwell.core.callcaching._
import cromwell.core.io._
import cromwell.core.logging.JobLogging
import wdl4s.values.WdlFile

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
}

/** A default implementation of the cache hit copying params. */
case class DefaultStandardFileHashingActorParams
(
  override val jobDescriptor: BackendJobDescriptor,
  override val backendInitializationDataOption: Option[BackendInitializationData],
  override val serviceRegistryActor: ActorRef,
  override val ioActor: ActorRef,
  override val configurationDescriptor: BackendConfigurationDescriptor
) extends StandardFileHashingActorParams

class DefaultStandardFileHashingActor(standardParams: StandardFileHashingActorParams) extends StandardFileHashingActor(standardParams) with DefaultIoCommandBuilder

object StandardFileHashingActor {
  case class FileHashingFunction(work: (SingleFileHashRequest, LoggingAdapter) => Try[String])

  sealed trait BackendSpecificHasherCommand { def jobKey: JobKey }
  final case class SingleFileHashRequest(jobKey: JobKey, hashKey: HashKey, file: WdlFile, initializationData: Option[BackendInitializationData]) extends BackendSpecificHasherCommand

  sealed trait BackendSpecificHasherResponse extends SuccessfulHashResultMessage
  case class FileHashResponse(hashResult: HashResult) extends BackendSpecificHasherResponse { override def hashes = Set(hashResult) }
}

abstract class StandardFileHashingActor(standardParams: StandardFileHashingActorParams) extends Actor with ActorLogging with JobLogging with IoClientHelper with StandardCachingActorHelper {
  this: IoCommandBuilder =>
  override val ioActor = standardParams.ioActor
  override lazy val jobDescriptor: BackendJobDescriptor = standardParams.jobDescriptor
  override lazy val backendInitializationDataOption: Option[BackendInitializationData] = standardParams.backendInitializationDataOption
  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor
  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor

  def customHashStrategy(fileRequest: SingleFileHashRequest): Option[Try[String]] = None

  def fileHashingReceive: Receive = {
    // Hash Request
    case fileRequest: SingleFileHashRequest =>
      val replyTo = sender()

      customHashStrategy(fileRequest) match {
        case Some(Success(result)) => context.parent ! FileHashResponse(HashResult(fileRequest.hashKey, HashValue(result)))
        case Some(Failure(failure)) => context.parent ! HashingFailedMessage(fileRequest.file.value, failure)
        case None => asyncHashing(fileRequest, replyTo)
      }

    // Hash Success
    case (fileHashRequest: SingleFileHashRequest, IoSuccess(_, result: String)) =>
      context.parent ! FileHashResponse(HashResult(fileHashRequest.hashKey, HashValue(result)))

    // Hash Failure
    case (fileHashRequest: SingleFileHashRequest, IoFailure(_, failure: Throwable)) =>
      context.parent ! HashingFailedMessage(fileHashRequest.file.value, failure)

    case other =>
      log.warning(s"Async File hashing actor received unexpected message: $other")
  }

  def asyncHashing(fileRequest: SingleFileHashRequest, replyTo: ActorRef) = getPath(fileRequest.file.value) match {
    case Success(gcsPath) => sendIoCommandWithContext(hashCommand(gcsPath), fileRequest)
    case Failure(failure) => replyTo ! HashingFailedMessage(fileRequest.file.value, failure)
  }

  override def receive: Receive = ioReceive orElse fileHashingReceive

  protected def onTimeout(message: Any, to: ActorRef): Unit = {
    context.parent ! HashingServiceUnvailable
  }
}
