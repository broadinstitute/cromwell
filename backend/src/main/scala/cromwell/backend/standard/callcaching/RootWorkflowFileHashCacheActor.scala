package cromwell.backend.standard.callcaching

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import com.google.common.cache.{CacheBuilder, CacheLoader, LoadingCache}
import cromwell.backend.standard.callcaching.RootWorkflowFileHashCacheActor.IoHashCommandWithContext
import cromwell.core.actor.RobustClientHelper.RequestTimeout
import cromwell.core.io._


class RootWorkflowFileHashCacheActor private(override val ioActor: ActorRef) extends Actor with ActorLogging with IoClientHelper {
  case class FileHashRequester(replyTo: ActorRef, fileHashContext: FileHashContext, ioCommand: IoCommand[_])

  sealed trait FileHashValue
  // The hash value is not yet in the cache and has not been requested.
  case object FileHashValueNotRequested extends FileHashValue
  // The hash value has been requested but is not yet in the cache.
  case class FileHashValueRequested(requesters: NonEmptyList[FileHashRequester]) extends FileHashValue
  // Hashing succeeded.
  case class FileHashSuccess(value: String) extends FileHashValue
  // Hashing failed.
  case class FileHashFailure(error: String) extends FileHashValue

  val cache: LoadingCache[String, FileHashValue] = CacheBuilder.newBuilder().build(
    new CacheLoader[String, FileHashValue] {
      override def load(key: String): FileHashValue = FileHashValueNotRequested
    })

  override def receive: Receive = ioReceive orElse cacheOrHashReceive

  val cacheOrHashReceive: Receive = {
    // Hash Request
    case hashCommand: IoHashCommandWithContext =>
      val key = hashCommand.fileHashContext.file
      lazy val requester = FileHashRequester(sender, hashCommand.fileHashContext, hashCommand.ioHashCommand)
      cache.get(key) match {
        case FileHashValueNotRequested =>
          // The hash is not in the cache and has not been requested. Make the hash request and register this requester
          // to be notified when the hash value becomes available.
          sendIoCommandWithContext(hashCommand.ioHashCommand, hashCommand.fileHashContext)
          cache.put(key, FileHashValueRequested(requesters = NonEmptyList.of(requester)))
        case FileHashValueRequested(requesters) =>
          // We don't have the hash but it has already been requested. Just add this requester and continue waiting for the
          // hash to become available.
          cache.put(key, FileHashValueRequested(requesters = requester :: requesters))
        case FileHashSuccess(value) =>
          sender ! Tuple2(hashCommand.fileHashContext, IoSuccess(requester.ioCommand, value))
        case FileHashFailure(error) =>
          sender ! Tuple2(hashCommand.fileHashContext, IoFailure(requester.ioCommand, new IOException(error)))
      }
    // Hash Success
    case (hashContext: FileHashContext, success @ IoSuccess(_, value: String)) =>
      handleHashResult(success, hashContext) { requesters =>
        requesters.toList foreach { case FileHashRequester(replyTo, fileHashContext, ioCommand) =>
          replyTo ! Tuple2(fileHashContext, IoSuccess(ioCommand, success.result))
        }
        cache.put(hashContext.file, FileHashSuccess(value))
      }
    // Hash Failure
    case (hashContext: FileHashContext, failure: IoFailure[_]) =>
      handleHashResult(failure, hashContext) { requesters =>
        requesters.toList foreach { case FileHashRequester(replyTo, fileHashContext, ioCommand) =>
          replyTo ! Tuple2(fileHashContext, IoFailure(ioCommand, failure.failure))
        }
        cache.put(hashContext.file, FileHashFailure(s"Error hashing file '${hashContext.file}': ${failure.failure.getMessage}"))
      }
    case other =>
      log.warning(s"Root workflow file hash caching actor received unexpected message: $other")
  }

  // Invoke the supplied block on the happy path, handle unexpected states for IoSuccess and IoFailure with common code.
  private def handleHashResult(ioAck: IoAck[_], fileHashContext: FileHashContext)
                              (notifyRequestersAndCacheValue: NonEmptyList[FileHashRequester] => Unit): Unit = {
    cache.get(fileHashContext.file) match {
      case FileHashValueRequested(requesters) => notifyRequestersAndCacheValue(requesters)
      case FileHashValueNotRequested =>
        log.error(s"Programmer error! Not expecting message type ${ioAck.getClass.getSimpleName} with no requesters for the hash: $fileHashContext")
      case _ =>
        log.error(s"Programmer error! Not expecting message type ${ioAck.getClass.getSimpleName} when the hash value has already been received: $fileHashContext")
    }
  }

  override protected def onTimeout(message: Any, to: ActorRef): Unit = {
    message match {
      case (fileHashContext: FileHashContext, _) =>
        // Send this message to all requestors.
        cache.get(fileHashContext.file) match {
          case FileHashValueRequested(requesters) =>
            requesters.toList foreach { case FileHashRequester(replyTo, requestContext, ioCommand) => replyTo ! RequestTimeout(Tuple2(requestContext, ioCommand), replyTo) }
            // Allow for the possibility of trying again on a timeout.
            cache.put(fileHashContext.file, FileHashValueNotRequested)
          case FileHashValueNotRequested =>
            log.error(s"Programmer error! Not expecting a hash request timeout when a hash value has not been requested: ${fileHashContext.file}")
          case _ =>
            log.error(s"Programmer error! Not expecting a hash request timeout when a hash value is already in the cache: ${fileHashContext.file}")
        }
      case other =>
        log.error(s"Programmer error! Root workflow file hash caching actor received unexpected timeout message: $other")
    }
  }
}

object RootWorkflowFileHashCacheActor {
  case class IoHashCommandWithContext(ioHashCommand: IoHashCommand, fileHashContext: FileHashContext)

  def props(ioActor: ActorRef): Props = Props(new RootWorkflowFileHashCacheActor(ioActor))
}
