package cromwell.database.sql

import java.sql.{Blob, Clob, Timestamp}
import java.time.{OffsetDateTime, ZoneId}
import javax.sql.rowset.serial.{SerialBlob, SerialClob}

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

  implicit class ClobOptionToRawString(val clobOption: Option[Clob]) extends AnyVal {
    // yes, it starts at 1
    def toRawStringOption: Option[String] = clobOption.map(clob => clob.getSubString(1, clob.length.toInt))

    def toRawString: String = toRawStringOption.getOrElse("")

    def parseSystemTimestampOption: Option[Timestamp] = toRawStringOption map { rawString =>
      OffsetDateTime.parse(rawString).toSystemTimestamp
    }
  }

  implicit class StringOptionToClobOption(val strOption: Option[String]) extends AnyVal {
    def toClob: Option[Clob] = strOption.flatMap(_.toClob)
  }

  implicit class StringToClobOption(val str: String) extends AnyVal {
    def toClob: Option[Clob] = if (str.isEmpty) None else Option(new SerialClob(str.toCharArray))
  }

  implicit class BlobToBytes(val blobOption: Option[Blob]) extends AnyVal {
    // yes, it starts at 1
    def toBytesOption: Option[Array[Byte]] = blobOption.map(blob => blob.getBytes(1, blob.length.toInt))

    def toBytes: Array[Byte] = toBytesOption.getOrElse(Array.empty)
  }

  implicit class BytesOptionToBlob(val bytesOption: Option[Array[Byte]]) extends AnyVal {
    def toBlob: Option[Blob] = bytesOption.flatMap(_.toBlob)
  }

  implicit class BytesToBlob(val bytes: Array[Byte]) extends AnyVal {
    def toBlob: Option[Blob] = if (bytes.isEmpty) None else Option(new SerialBlob(bytes))
  }
}
