package cromwell.core.io

import akka.actor.{Actor, ActorLogging, ActorRef}
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path

import scala.concurrent.duration._
import scala.concurrent.{Future, Promise}
import scala.language.postfixOps


trait AsyncIo extends IoClientHelper { this: Actor with ActorLogging with IoCommandBuilder =>
  
  protected val ioTimeout = 3 minutes
  
  private final def ioResponseReceive: Receive = {
    case (promise: Promise[_], ack: IoAck[Any] @ unchecked) =>
      cancelTimeout(promise -> ack.command)
      // This is not typesafe. 
      // However the sendIoCommand method ensures that the command and the promise have the same generic type
      // Which means as long as only the sendIoCommand method is used to send requests, and the ioActor honors his contract
      // and send back the right context with the right response, the types are virtually guaranteed to match.
      promise.asInstanceOf[Promise[Any]].complete(ack.toTry)
      ()
  }
  
  protected final def ioReceive = ioResponseReceive orElse robustReceive

  /**
    * IMPORTANT: This loads the entire content of the file into memory !
    * Only use for small files !
    */
  def contentAsStringAsync(path: Path): Future[String] = {
    val promise = Promise[String]
    sendIoCommandWithPromise(contentAsStringCommand(path), promise)
    promise.future
  }
  
  def writeAsync(path: Path, content: String, options: OpenOptions): Future[Unit] = {
    val promise = Promise[Unit]
    sendIoCommandWithPromise(writeCommand(path, content, options), promise)
    promise.future
  }

  def sizeAsync(path: Path): Future[Long] = {
    val promise = Promise[Long]
    sendIoCommandWithPromise(sizeCommand(path), promise)
    promise.future
  }

  def deleteAsync(path: Path, swallowIoExceptions: Boolean = false): Future[Unit] = {
    val promise = Promise[Unit]
    sendIoCommandWithPromise(deleteCommand(path, swallowIoExceptions), promise)
    promise.future
  }

  def copyAsync(src: Path, dest: Path, overwrite: Boolean = true): Future[Unit] = {
    val promise = Promise[Unit]
    // Allow for a much larger timeout for copies, as large files can take a while (even on gcs, if they are in different locations...)
    sendIoCommandWithPromise(copyCommand(src, dest, overwrite), promise, 1 hour)
    promise.future
  }

  private def sendIoCommandWithPromise[T](command: IoCommand[T], promise: Promise[T], timeout: FiniteDuration = ioTimeout) = {
    sendIoCommandWithContext(command, promise, timeout)
  }

  override def onTimeout(message: Any, to: ActorRef): Unit = message match {
    case (promise: Promise[_], ioAck: IoAck[_]) => 
      promise.tryFailure(IoTimeout(ioAck.command))
      ()
    case _ =>
  }
}
