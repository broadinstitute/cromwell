package cromwell.core.logging

import ch.qos.logback.classic.pattern.ThreadConverter
import ch.qos.logback.classic.spi.ILoggingEvent

/**
  * Log the Akka sourceThread if found, otherwise log the event thread.
  * 
  *   - https://doc.akka.io/docs/akka/current/logging.html#logging-thread-akka-source-and-actor-system-in-mdc
  *   - https://logback.qos.ch/manual/layouts.html#customConversionSpecifier
  */
class EnhancedThreadConverter extends ThreadConverter {
  override def convert(event: ILoggingEvent): String = {
    val mdc = event.getMDCPropertyMap
    if (mdc.containsKey("sourceThread")) {
      mdc.get("sourceThread")
    } else {
      super.convert(event)
    }
  }
}
