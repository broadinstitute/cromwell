package cromwell.database.sql

import java.sql.{Clob, Timestamp}
import java.time.{OffsetDateTime, ZoneId}
import javax.sql.rowset.serial.SerialClob

object SqlConverters {

  // TODO: Storing times relative to system zone. Look into db/slick using OffsetDateTime, or storing datetimes as UTC?
  // http://stackoverflow.com/questions/34608650/scala-slick-3-0-implicit-mapping-between-java8-offsetdatetime-and-timestamp
  // https://github.com/slick/slick/issues/1026

  implicit class TimestampToSystemOffsetDateTime(val timestamp: Timestamp) extends AnyVal {
    def toSystemOffsetDateTime = timestamp.toLocalDateTime.atZone(ZoneId.systemDefault).toOffsetDateTime
  }

  implicit class OffsetDateTimeToSystemTimestamp(val offsetDateTime: OffsetDateTime) extends AnyVal {
    def toSystemTimestamp = Timestamp.valueOf(offsetDateTime.atZoneSameInstant(ZoneId.systemDefault).toLocalDateTime)
  }

  implicit class ClobToRawString(val clob: Clob) extends AnyVal {
    def toRawString: String = clob.getSubString(1, clob.length.toInt) // yes, it starts at 1
  }

  implicit class StringToClob(val str: String) extends AnyVal {
    def toClob: Clob = new SerialClob(str.toCharArray)

    def toNonEmptyClob: Option[Clob] = if (str.isEmpty) None else Option(new SerialClob(str.toCharArray))
  }

}
