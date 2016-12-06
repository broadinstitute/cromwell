package cromwell.logging

import java.time.OffsetDateTime
import java.time.format.DateTimeFormatter

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.pattern.ThrowableProxyConverter
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.LayoutBase
import lenthall.util.TerminalUtil

object TerminalLayout {
  val Converter = new ThrowableProxyConverter
  Converter.start()

  private val dateTimeFormatter = DateTimeFormatter.ofPattern("uuuu-MM-dd HH:mm:ss,SS")

  implicit class EnhancedILoggingEvent(val event: ILoggingEvent) extends AnyVal {
    def toStackTrace: String = Converter.convert(event)
  }

  implicit class ColorString(msg: String) {
    def colorizeUuids: String = {
      "UUID\\((.*?)\\)".r.findAllMatchIn(msg).foldLeft(msg) {
        case (l, r) =>
          val color = if (Option(System.getProperty("RAINBOW_UUID")).isDefined)
            Math.abs(17 * r.group(1).substring(0,8).map(_.toInt).product) % 209 + 22
          else 2
          l.replace(r.group(0), TerminalUtil.highlight(color, r.group(1)))
      }
    }
    def colorizeCommand: String = msg.replaceAll("`([^`]*?)`", TerminalUtil.highlight(5, "$1"))
  }
}

class TerminalLayout extends LayoutBase[ILoggingEvent] {
  import TerminalLayout._
  def doLayout(event: ILoggingEvent): String = {
    val level = event.getLevel match {
      case Level.WARN => TerminalUtil.highlight(220, "warn")
      case Level.ERROR => TerminalUtil.highlight(1, "error")
      case x => x.toString.toLowerCase
    }
    val timestamp = OffsetDateTime.now.format(dateTimeFormatter)
    val highlightedMessage = event.getFormattedMessage.colorizeUuids.colorizeCommand
    s"[$timestamp] [$level] $highlightedMessage\n${event.toStackTrace}"
  }
}
