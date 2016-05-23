package cromwell.core

import java.sql.Timestamp
import java.time.{ZoneId, OffsetDateTime}
import java.time.format._
import java.time.format._

object KnowsWhatTimeItIs {
  implicit class JodaTimestampFormatter(val underlying: Timestamp) extends AnyVal {
    def asJodaString = underlying.toLocalDateTime.atZone(ZoneId.systemDefault).format(DateTimeFormatter.ISO_OFFSET_DATE_TIME)
    def toOffsetDateTime = underlying.toLocalDateTime.atZone(ZoneId.systemDefault).toOffsetDateTime
  }
}

trait KnowsWhatTimeItIs {
  def currentTime = Timestamp.valueOf(OffsetDateTime.now.toLocalDateTime)
}
