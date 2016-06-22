package cromwell.logging

import java.time.OffsetDateTime
import java.time.format.DateTimeFormatter

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.pattern.ThrowableProxyConverter
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.LayoutBase

object TerminalLayout {
  val HighlightedError = fansi.Color.Red("error")
  val HighlightedWarn = fansi.Color.Yellow("warn")

  val Converter = new ThrowableProxyConverter
  Converter.start()

  private val dateTimeFormatter = DateTimeFormatter.ofPattern("uuuu-MM-dd HH:mm:ss,SS")

  implicit class EnhancedILoggingEvent(val event: ILoggingEvent) extends AnyVal {
    def toStackTrace: String = Converter.convert(event)
  }

  implicit class ColorString(msg: String) {
    def colorizeUuids: String = {
      "UUID\\((.*?)\\)".r.findAllMatchIn(msg).foldLeft(msg) {
        case (l, r) => l.replace(r.group(0), fansi.Color.Green(r.group(1)).render)
      }
    }

    def colorizeCommand: String = msg.replaceAll("`([^`]*?)`", fansi.Color.Magenta("$1").render)
  }
}

class TerminalLayout extends LayoutBase[ILoggingEvent] {
  import TerminalLayout._
  def doLayout(event: ILoggingEvent): String = {
    val level = event.getLevel match {
      case Level.WARN => HighlightedWarn
      case Level.ERROR => HighlightedError
      case x => x.toString.toLowerCase
    }

    val timestamp = OffsetDateTime.now.format(dateTimeFormatter)
    val highlightedMessage = event.getFormattedMessage.colorizeUuids.colorizeCommand
    s"[$timestamp] [$level] $highlightedMessage\n${event.toStackTrace}"
  }
}
