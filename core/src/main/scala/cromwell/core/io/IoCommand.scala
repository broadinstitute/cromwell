package cromwell.core.io

import java.time.OffsetDateTime
import java.util.UUID

import better.files.File.OpenOptions
import com.google.api.client.util.ExponentialBackOff
import common.util.Backoff
import common.util.StringUtil.EnhancedToStringable
import cromwell.core.io.IoContentAsStringCommand.IoReadOptions
import cromwell.core.path.Path
import cromwell.core.retry.SimpleExponentialBackoff
import org.slf4j.LoggerFactory

import scala.concurrent.duration._
import scala.language.postfixOps

object IoCommand {
  private val logger = LoggerFactory.getLogger(IoCommand.getClass)
  val IOCommandWarnLimit: FiniteDuration = 5 minutes

  def defaultGoogleBackoff: ExponentialBackOff = new ExponentialBackOff.Builder()
    .setInitialIntervalMillis((1 second).toMillis.toInt)
    .setMaxIntervalMillis((5 minutes).toMillis.toInt)
    .setMultiplier(3L)
    .setRandomizationFactor(0.2D)
    .setMaxElapsedTimeMillis((10 minutes).toMillis.toInt)
    .build()
  
  def defaultBackoff: Backoff = SimpleExponentialBackoff(defaultGoogleBackoff)
  
  type RetryCommand[T] = (FiniteDuration, IoCommand[T])
}

trait IoCommand[+T] {
  private val uuid = UUID.randomUUID().toString
  private val creation = OffsetDateTime.now

  def commandDescription: String

  def logIOMsgOverLimit(message: => String): Unit = {
    val millis: Long = java.time.Duration.between(creation, OffsetDateTime.now).toMillis
    if (millis > IoCommand.IOCommandWarnLimit.toMillis) {
      val seconds = millis / 1000D

      /*
        For now we decided to log this as INFO. In future if needed, we can update this to WARN.
        Note: Below links would be useful to look at in case we want to increase the visibility of this message
              (send an INFO message to Sentry) or even silence it all together (not even send it to Stackdriver):
              - Cromwell server `logback.xml` with various environment-variable based conditionals
                (https://github.com/broadinstitute/cromwell/blob/53/server/src/main/resources/logback.xml)
              - Environment variables being set to affect the server `logback.xml`
                (https://github.com/broadinstitute/firecloud-develop/blob/c77e0f371be0aac545e204f1a134cc6f8ef3c301/run-context/live/configs/cromwell/app.env.ctmpl#L42-L51)
              - Logback manual (http://logback.qos.ch/manual/index.html)
       */
      IoCommand.logger.info(f"(IO-$uuid) '$message' is over 5 minutes. It was running for " +
        f"$seconds%,.3f seconds. IO command description: '$commandDescription'")
    }
  }

  /**
    * Completes the command successfully
    * @return a message to be sent back to the sender, if needed
    */
  def success[S >: T](value: S): IoSuccess[S] = {
    logIOMsgOverLimit(s"IOCommand.success '$value'")
    IoSuccess(this, value)
  }
  
  /**
    * Fail the command with an exception
    */
  def fail[S >: T](failure: Throwable): IoFailure[S] = {
    logIOMsgOverLimit(s"IOCommand.fail '${failure.toPrettyElidedString(limit = 1000)}'")
    IoFailure(this, failure)
  }

  def failReadForbidden[S >: T](failure: Throwable, forbiddenPath: String): IoReadForbiddenFailure[S] = {
    logIOMsgOverLimit(s"IOCommand.failReadForbidden '${failure.toPrettyElidedString(limit = 1000)}' path '$forbiddenPath'")
    IoReadForbiddenFailure(this, failure, forbiddenPath)
  }

  /**
    * A short name describing the command
    */
  def name: String
}

trait SingleFileIoCommand[T] extends IoCommand[T] {
  def file: Path
}

/**
  * Copy source -> destination
  * Will create the destination directory if it doesn't exist.
  */
abstract class IoCopyCommand(val source: Path, val destination: Path) extends IoCommand[Unit] {
  override def toString = s"copy ${source.pathAsString} to ${destination.pathAsString}"
  override lazy val name = "copy"
}
  
object IoContentAsStringCommand {

  /**
    * Options to customize reading of a file.
    * @param maxBytes If specified, only reads up to maxBytes Bytes from the file
    * @param failOnOverflow If this is true, maxBytes is specified, and the file is larger than maxBytes, fail the command. 
    */
  case class IoReadOptions(maxBytes: Option[Int], failOnOverflow: Boolean)
}

/**
  * Read file as a string (load the entire content in memory)
  */
abstract class IoContentAsStringCommand(val file: Path, val options: IoReadOptions = IoReadOptions(None, failOnOverflow = false)) extends SingleFileIoCommand[String] {
  override def toString = s"read content of ${file.pathAsString}"
  override lazy val name = "read"
}

/**
  * Return the size of file
  */
abstract class IoSizeCommand(val file: Path) extends SingleFileIoCommand[Long] {
  override def toString = s"get size of ${file.pathAsString}"
  override lazy val name = "size"
}

/**
  * Write content in file
  * Will create the destination directory if it doesn't exist.
  */
abstract class IoWriteCommand(val file: Path,
                              val content: String,
                              val openOptions: OpenOptions,
                              val compressPayload: Boolean) extends SingleFileIoCommand[Unit] {
  override def toString = s"write to ${file.pathAsString}"
  override lazy val name = "write"
}

/**
  * Delete file
  */
abstract class IoDeleteCommand(val file: Path, val swallowIOExceptions: Boolean) extends SingleFileIoCommand[Unit] {
  override def toString = s"delete ${file.pathAsString}"
  override lazy val name = "delete"
}

/**
  * Get Hash value for file
  */
abstract class IoHashCommand(val file: Path) extends SingleFileIoCommand[String] {
  override def toString = s"get hash of ${file.pathAsString}"
  override lazy val name = "hash"
}

/**
  * Touch a file
  */
abstract class IoTouchCommand(val file: Path) extends SingleFileIoCommand[Unit] {
  override def toString = s"touch ${file.pathAsString}"
  override lazy val name = "touch"
}

/**
  * Check whether a file exists
  */
abstract class IoExistsCommand(val file: Path) extends SingleFileIoCommand[Boolean] {
  override def toString = s"check whether ${file.pathAsString} exists"
  override lazy val name = "exist"
}

/**
  * Return the lines of a file in a collection
  */
abstract class IoReadLinesCommand(val file: Path) extends SingleFileIoCommand[Traversable[String]] {
  override def toString = s"read lines of ${file.pathAsString}"
  override lazy val name = "read lines"
}

/**
  * Check whether a path represents a directory
  */
abstract class IoIsDirectoryCommand(val file: Path) extends SingleFileIoCommand[Boolean] {
  override def toString = s"check whether ${file.pathAsString} is a directory"
  override lazy val name = "is directory"
}
