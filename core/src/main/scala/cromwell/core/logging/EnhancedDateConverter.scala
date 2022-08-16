package cromwell.core.logging

import ch.qos.logback.classic.pattern.DateConverter
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.CoreConstants
import ch.qos.logback.core.util.CachingDateFormatter

import java.util.TimeZone
import scala.jdk.CollectionConverters._

/**
  * Log the Akka akkaTimestamp if found in the MDC, otherwise log the original event timestamp.
  *
  *   - https://doc.akka.io/docs/akka/current/logging.html#more-accurate-timestamps-for-log-output-in-mdc
  *   - https://logback.qos.ch/manual/layouts.html#customConversionSpecifier
  *
  * NOTE: For proper configuration both this EnhancedDateConverter should be configured into the logback.xml AND the
  * configuration file should set akka.loggers = ["cromwell.core.logging.EnhancedSlf4jLogger"].
  */
class EnhancedDateConverter extends DateConverter {
  protected var cachingDateFormatterProtected: CachingDateFormatter = _

  /* Duplicated from ch.qos.logback.classic.pattern.DateConverter as cachingDateFormatter is package private. */
  override def start(): Unit = {
    cachingDateFormatterProtected = Option(getFirstOption) match {
      case Some(CoreConstants.ISO8601_STR) | None => new CachingDateFormatter(CoreConstants.ISO8601_PATTERN)
      case Some(datePattern) =>
        try {
          new CachingDateFormatter(datePattern)
        } catch {
          case e: IllegalArgumentException =>
            addWarn("Could not instantiate SimpleDateFormat with pattern " + datePattern, e)
            // default to the ISO8601 format
            new CachingDateFormatter(CoreConstants.ISO8601_PATTERN)
        }
    }
    // if the option list contains a TZ option, then set it.
    Option(getOptionList)
      .toList
      .flatMap(_.asScala)
      .drop(1)
      .headOption
      .map(TimeZone.getTimeZone)
      .foreach(cachingDateFormatterProtected.setTimeZone)

    // Allow the parent class to start/initialize its private members.
    super.start()
  }

  /**
    * Look for the Akka timestamp and use that to format the date.
    *
    * Until this (currently 6+ year) issue is resolved, formatting the date as a Long requires using the
    * [[EnhancedSlf4jLogger]] versus Akka's basic Slf4jLogger.
    *
    *   - https://github.com/akka/akka/issues/18079#issuecomment-125175884
    */
  override def convert(event: ILoggingEvent): String = {
    val mdc = event.getMDCPropertyMap
    if (mdc.containsKey("akkaTimestamp")) {
      val timestamp = mdc.get("akkaTimestamp")
      timestamp.toLongOption match {
        case Some(value) => cachingDateFormatterProtected.format(value)
        case None => timestamp // Return the original timestamp string.
      }
    } else {
      super.convert(event)
    }
  }
}
