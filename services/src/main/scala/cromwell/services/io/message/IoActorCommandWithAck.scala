package cromwell.services.io.message

import java.nio.file.Path

import better.files.File.OpenOptions
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.services.io.IoActorCommand._
import cromwell.services.io._

sealed trait IoActorCommandWithAck[T] { this: IoActorCommand[T] =>
  override def success(value: T) = Option(IoSuccess(this, value))
  override def fail(failure: Throwable) = Option(IoFailure(this, failure))
  override def retried(failure: Throwable) = Option(IoRetried(this, failure))
}

/**
  * Copy source -> destination
  */
case class CopyCommandMessage(source: Path,
                              destination: Path,
                              overwrite: Boolean = true,
                              backoff: SimpleExponentialBackoff = defaultBackoff
                             ) extends CopyCommand with IoActorCommandWithAck[Unit] {
  override def withNextBackoff: CopyCommand = copy(backoff = backoff.next)
}

/**
  * Read file as a string (load the entire content in memory)
  */
case class ReadAsStringCommandMessage(file: Path,
                                      backoff: SimpleExponentialBackoff = defaultBackoff
                                     ) extends ReadAsStringCommand with IoActorCommandWithAck[String] {
  override def withNextBackoff: ReadAsStringCommand = copy(backoff = backoff.next)
}

/**
  * Return the size of file
  */
case class SizeCommandMessage(file: Path,
                              backoff: SimpleExponentialBackoff = defaultBackoff
                             ) extends SizeCommand with IoActorCommandWithAck[Long] {
  override def withNextBackoff = copy(backoff = backoff.next)
}

/**
  * Write content in file
  */
case class WriteCommandMessage(file: Path,
                               content: String,
                               openOptions: OpenOptions,
                               backoff: SimpleExponentialBackoff = defaultBackoff
                              ) extends WriteCommand with IoActorCommandWithAck[Unit] {
  override def withNextBackoff = copy(backoff = backoff.next)
}

/**
  * Delete file
  */
case class DeleteCommandMessage(file: Path,
                                swallowIOExceptions: Boolean = false,
                                backoff: SimpleExponentialBackoff = defaultBackoff
                               ) extends DeleteCommand with IoActorCommandWithAck[Unit] {
  override def withNextBackoff = copy(backoff = backoff.next)
}