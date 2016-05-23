package cromwell.core

import java.sql.Timestamp
import java.time.OffsetDateTime
import org.joda.time.DateTime

object KnowsWhatTimeItIs {
  implicit class JodaTimestampFormatter(val underlying: Timestamp) extends AnyVal {
    def asJodaString = new DateTime(underlying).toString
  }
}

trait KnowsWhatTimeItIs {
  def currentTime = Timestamp.valueOf(OffsetDateTime.now.toLocalDateTime)
}
