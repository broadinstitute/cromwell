package cromwell.engine.io

import com.typesafe.scalalogging.StrictLogging
import cromwell.core.io.IoCommand

import java.time.OffsetDateTime
import java.time.temporal.ChronoUnit
import scala.concurrent.duration.FiniteDuration

trait IoCommandStalenessBackpressuring extends StrictLogging {

  def maxStaleness: FiniteDuration

  private def commandStalenessThreshold: OffsetDateTime = OffsetDateTime.now().minusSeconds(maxStaleness.toSeconds)

  private def logAndBackpressure(ioCommand: IoCommand[_], onBackpressure: Option[Double] => Unit): Unit = {
    val millis = ChronoUnit.MILLIS.between(ioCommand.creation, commandStalenessThreshold)
    val scale = (maxStaleness.toMillis + millis.toDouble) / maxStaleness.toMillis

    val seconds = millis / 1000.0
    logger.info("I/O command {} seconds stale, applying I/O subsystem backpressure with scale {}",
      f"$seconds%,.3f", f"$scale%.2f")

    onBackpressure(Option(scale))
  }

  /** Invokes `onBackpressure` if `ioCommand` is older than the staleness limit returned by `maxStaleness`. */
  def backpressureIfStale(ioCommand: IoCommand[_], onBackpressure: Option[Double] => Unit): Unit = {
    if (ioCommand.creation.isBefore(commandStalenessThreshold)) {
      logAndBackpressure(ioCommand, onBackpressure)
    }
  }

  /** Invokes `onBackpressure` if at least one IoCommandContext in `contexts` is older than the
    * staleness limit returned by `maxStaleness`. */
  def backpressureIfStale(contexts: Seq[IoCommandContext[_]], onBackpressure: Option[Double] => Unit): Unit = {
    val oldest = contexts.minBy(_.creationTime)
    if (oldest.creationTime.isBefore(commandStalenessThreshold)) {
      logAndBackpressure(oldest.request, onBackpressure)
    }
  }
}
