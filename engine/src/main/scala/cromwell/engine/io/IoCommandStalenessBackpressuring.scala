package cromwell.engine.io

import com.typesafe.scalalogging.StrictLogging
import cromwell.core.io.IoCommand

import java.time.OffsetDateTime
import scala.concurrent.duration.FiniteDuration

trait IoCommandStalenessBackpressuring extends StrictLogging {

  def maxStaleness: FiniteDuration

  private def commandStalenessThreshold: OffsetDateTime = OffsetDateTime.now().minusSeconds(maxStaleness.toSeconds)

  /** Invokes the specified `onBackpressure` function if the specified `ioCommand` is older than the staleness
    * limit returned by `maxStaleness`. */
  def backpressureIfStale(ioCommand: IoCommand[_], onBackpressure: () => Unit): Unit = {
    if (ioCommand.creation.isBefore(commandStalenessThreshold)) {
      logger.info("Command with creation time {} staler than {}, applying I/O subsystem backpressure",
        ioCommand.creation, maxStaleness)
      onBackpressure()
    }
  }

  /** Invokes the specified `onBackpressure` function if at least one IoCommandContext in `contexts` is older than the
    * staleness limit returned by `maxStaleness`. */
  def backpressureIfStale(contexts: Seq[IoCommandContext[_]], onBackpressure: () => Unit): Unit = {
    val threshold = commandStalenessThreshold
    if (contexts.exists(_.creationTime.isBefore(threshold))) {
      logger.info("Found at least one I/O command staler than {}, applying I/O subsystem backpressure",
        maxStaleness)
      onBackpressure()
    }
  }
}
