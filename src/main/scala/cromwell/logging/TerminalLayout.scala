package cromwell.logging

import java.text.SimpleDateFormat
import java.util.Date

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.pattern.ThrowableProxyConverter
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.LayoutBase
import cromwell.logging.TerminalLayout.EnhancedILoggingEvent
import cromwell.util.TerminalUtil

class TerminalLayout extends LayoutBase[ILoggingEvent] {
  def doLayout(event: ILoggingEvent): String = {
    val level = event.getLevel match {
      case Level.WARN => TerminalUtil.highlight(220, "warn")
      case Level.ERROR => TerminalUtil.highlight(1, "error")
      case x => x.toString.toLowerCase
    }

    val timestamp = new SimpleDateFormat("YYYY-MM-dd HH:mm:ss,SS").format(new Date())

    val ansiColor = Option(System.getProperty("RAINBOW_UUID")) match {
      case Some(_) => "UUID\\((.*?)\\)".r.findFirstMatchIn(event.getFormattedMessage).map(s => (Math.abs(17 * s.group(1).map(_.toInt).product) % 209) + 22).getOrElse(2)
      case None => 2
    }

    val highlightedMessage = event.getFormattedMessage
      .replaceAll("UUID\\((.*?)\\)", TerminalUtil.highlight(ansiColor, "$1"))
      .replaceAll("`([^`]*?)`", TerminalUtil.highlight(5, "$1"))

    s"[$timestamp] [$level] $highlightedMessage\n${event.toStackTrace}"
  }
}

object TerminalLayout {
  val Converter = new ThrowableProxyConverter
  Converter.start()

  implicit class EnhancedILoggingEvent(val event: ILoggingEvent) extends AnyVal {
    def toStackTrace: String = Converter.convert(event)
  }
}
