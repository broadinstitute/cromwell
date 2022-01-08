package common.util

import java.time.format.DateTimeFormatter
import java.time.{OffsetDateTime, ZoneOffset}

object TimeUtil {
  /**
    * Instead of "one of" the valid ISO-8601 formats, standardize on this one:
    * https://github.com/openjdk/jdk/blob/jdk8-b120/jdk/src/share/classes/java/time/OffsetDateTime.java#L1886
    */
  private val Iso8601MillisecondsFormat = DateTimeFormatter.ofPattern("uuuu-MM-dd'T'HH:mm:ss.SSSXXXXX")

  implicit class EnhancedOffsetDateTime(val offsetDateTime: OffsetDateTime) extends AnyVal {
    /**
      * Discards the original timezone and shifts the time to UTC, then returns the ISO-8601 formatted string with
      * exactly three digits of milliseconds.
      */
    def toUtcMilliString: String = Option(offsetDateTime).map(
      _.atZoneSameInstant(ZoneOffset.UTC).format(Iso8601MillisecondsFormat)
    ).orNull
  }
}
