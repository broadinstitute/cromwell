package cromwell.services.io

import java.nio.file.Path

import better.files.File.OpenOptions
import com.google.api.client.util.{BackOff, ExponentialBackOff}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

import scala.concurrent.duration.{FiniteDuration, _}
import scala.language.postfixOps

object IoActorCommand {
  def defaultGoogleBackoff = new ExponentialBackOff.Builder()
    .setInitialIntervalMillis((1 second).toMillis.toInt)
    .setMaxIntervalMillis((5 minutes).toMillis.toInt)
    .setMultiplier(5L)
    .setMaxElapsedTimeMillis((10 minutes).toMillis.toInt)
    .build()
  def defaultBackoff = SimpleExponentialBackoff(defaultGoogleBackoff)
}

/**
  * Io Command
  *
  * @tparam T type of the value to be returned after the command is executed
  */
sealed trait IoActorCommand[T] extends ServiceRegistryMessage {
  override val serviceName = IoService.serviceName

  def backoff: SimpleExponentialBackoff

  /**
    * Current value of the backoff duration or None if the command should not be retried
    */
  def timeBeforeNextTry: Option[FiniteDuration] = backoff.backoffMillis match {
    case BackOff.STOP => None
    case value => Option(value millis)
  }

  /**
    * Increment the exponential backoff value
    */
  def withNextBackoff: IoActorCommand[T]

  /**
    * Fail the command with an exception
    */
  def fail(failure: Throwable): Option[Any]

  /**
    * Message emitted when the command failed but is being re-submitted for retry
    */
  def retried(failure: Throwable): Option[Any]

  /**
    * Completes the command successfully
    * @return a message to be sent back to the sender, if needed
    */
  def success(value: T): Option[Any]
}

/**
  * Copy source -> destination
  */
trait CopyCommand extends IoActorCommand[Unit] {
  def source: Path
  def destination: Path
  def overwrite: Boolean
} 

/**
  * Read file as a string (load the entire content in memory)
  */
trait ReadAsStringCommand extends IoActorCommand[String] {
  def file: Path
} 

/**
  * Return the size of file
  */
trait SizeCommand extends IoActorCommand[Long] {
  def file: Path
}

/**
  * Write content in file
  */
trait WriteCommand extends IoActorCommand[Unit] {
  def file: Path
  def content: String
  def openOptions: OpenOptions
}

/**
  * Delete file
  */
trait DeleteCommand extends IoActorCommand[Unit] {
  def file: Path
  def swallowIOExceptions: Boolean
}
