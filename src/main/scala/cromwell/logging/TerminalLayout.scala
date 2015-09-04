package cromwell.logging

import java.util

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.pattern.ThrowableProxyConverter
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.{ConsoleAppender, LayoutBase}
import cromwell.util.TerminalUtil
import TerminalLayout.EnhancedILoggingEvent

class TerminalLayout extends LayoutBase[ILoggingEvent] {
  def doLayout(event: ILoggingEvent): String = {
    val level = event.getLevel match {
      case Level.WARN => TerminalUtil.highlight(220, "warn")
      case Level.ERROR => TerminalUtil.highlight(1, "error")
      case x => x.toString.toLowerCase
    }

    val highlightedMessage = event.getFormattedMessage
      .replaceAll("UUID\\((.*?)\\)", TerminalUtil.highlight(2, "$1"))
      .replaceAll("`([^`]*?)`", TerminalUtil.highlight(5, "$1"))

    /* For some reason a '{}' is the value of getMessage only for Actors.
       This prepends a highlighted asterisk to messages that come from actors.
     */
    val prefix = if (event.getMessage == "{}") s"[${TerminalUtil.highlight(129, "*")}] " else ""

    s"$prefix[$level] $highlightedMessage\n${event.toStackTrace}"
  }
}

object TerminalLayout {
  val Converter = new ThrowableProxyConverter
  Converter.start()

  implicit class EnhancedILoggingEvent(val event: ILoggingEvent) extends AnyVal {
    def toStackTrace: String = Converter.convert(event)
  }
}
