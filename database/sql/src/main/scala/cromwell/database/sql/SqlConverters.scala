package cromwell.database.sql

import java.sql.{Blob, Clob, Timestamp}
import java.time.{OffsetDateTime, ZoneId}
import javax.sql.rowset.serial.{SerialBlob, SerialClob}

import eu.timepit.refined.api.Refined
import eu.timepit.refined.collection.NonEmpty

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
    def toRawStringOption: Option[String] = clobOption.map(_.toRawString)

    def toRawString: String = toRawStringOption.getOrElse("")

    def parseSystemTimestampOption: Option[Timestamp] = toRawStringOption map { rawString =>
      OffsetDateTime.parse(rawString).toSystemTimestamp
    }
  }

  implicit class ClobToRawString(val clob: Clob) extends AnyVal {
    // yes, it starts at 1
    def toRawString: String = clob.getSubString(1, clob.length.toInt)
  }

  implicit class StringOptionToClobOption(val strOption: Option[String]) extends AnyVal {
    def toClobOption: Option[Clob] = strOption.flatMap(_.toClobOption)
  }

  implicit class StringToClobOption(val str: String) extends AnyVal {
    def toClobOption: Option[Clob] = if (str.isEmpty) None else Option(new SerialClob(str.toCharArray))

    def toClob(default: String Refined NonEmpty): Clob = new SerialClob(default.toString.toCharArray)
  }

  implicit class BlobToBytes(val blob: Blob) extends AnyVal {
    // yes, it starts at 1
    def toBytes: Array[Byte] = blob.getBytes(1, blob.length.toInt)
  }

  implicit class BlobOptionToBytes(val blobOption: Option[Blob]) extends AnyVal {
    def toBytesOption: Option[Array[Byte]] = blobOption.map(_.toBytes)

    def toBytes: Array[Byte] = toBytesOption.getOrElse(Array.empty)
  }

  implicit class BytesOptionToBlob(val bytesOption: Option[Array[Byte]]) extends AnyVal {
    def toBlobOption: Option[Blob] = bytesOption.flatMap(_.toBlobOption)
  }

  implicit class BytesToBlob(val bytes: Array[Byte]) extends AnyVal {
    def toBlobOption: Option[Blob] = if (bytes.isEmpty) None else Option(new SerialBlob(bytes))
  }
}
