package cromwell.services.io.promise

import java.nio.file.Path

import better.files.File._
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.services.io.IoActorCommand._
import cromwell.services.io._

import scala.concurrent.Promise


sealed trait IoActorPromiseCommand[T] { this: IoActorCommand[T] =>
  def promise: Promise[T]
  override def fail(failure: Throwable) = {
    promise.tryFailure(failure)
    None
  }
  override def success(value: T) = {
    promise.trySuccess(value)
    None
  }

  override def retried(failure: scala.Throwable) = None
}

/**
  * Copy source -> destination
  */
case class CopyCommandPromise(source: Path,
                              destination: Path,
                              overwrite: Boolean = true,
                              backoff: SimpleExponentialBackoff = defaultBackoff)(val promise: Promise[Unit]) extends IoActorPromiseCommand[Unit] with CopyCommand {
  override def withNextBackoff: CopyCommandPromise = copy(backoff = backoff.next)(promise)
}

/**
  * Read file as a string (load the entire content in memory)
  */
case class ReadAsStringCommandPromise(file: Path,
                                      backoff: SimpleExponentialBackoff = defaultBackoff)(val promise: Promise[String]) extends IoActorPromiseCommand[String] with ReadAsStringCommand {
  override def withNextBackoff: ReadAsStringCommandPromise = copy(backoff = backoff.next)(promise)
}

/**
  * Return the size of file
  */
case class SizeCommandPromise(file: Path,
                              backoff: SimpleExponentialBackoff = defaultBackoff)(val promise: Promise[Long]) extends IoActorPromiseCommand[Long] with SizeCommand {
  override def withNextBackoff: SizeCommandPromise = copy(backoff = backoff.next)(promise)
}

/**
  * Write content in file
  */
case class WriteCommandPromise(file: Path,
                               content: String,
                               openOptions: OpenOptions,
                               backoff: SimpleExponentialBackoff = defaultBackoff)(val promise: Promise[Unit]) extends IoActorPromiseCommand[Unit] with WriteCommand {
  override def withNextBackoff: WriteCommandPromise = copy(backoff = backoff.next)(promise)
}

/**
  * Delete file
  */
case class DeleteCommandPromise(file: Path,
                                swallowIOExceptions: Boolean,
                                backoff: SimpleExponentialBackoff = defaultBackoff)(val promise: Promise[Unit]) extends IoActorPromiseCommand[Unit] with DeleteCommand {
  override def withNextBackoff: DeleteCommandPromise = copy(backoff = backoff.next)(promise)
}