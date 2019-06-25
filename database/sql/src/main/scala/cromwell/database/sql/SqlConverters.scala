package cromwell.database.sql

import java.sql.{Blob, Clob, Timestamp}
import java.time.{OffsetDateTime, ZoneId}

import javax.sql.rowset.serial.{SerialBlob, SerialClob}

import scala.concurrent.duration.FiniteDuration

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
    def toRawString: String = {
      // See notes on empty clob issues in StringToClobOption
      val length = clob.length.toInt
      // yes, it starts at 1
      if (length == 0) "" else clob.getSubString(1, length)
    }
  }

  implicit class StringOptionToClobOption(val strOption: Option[String]) extends AnyVal {
    def toClobOption: Option[SerialClob] = strOption.flatMap(_.toClobOption)
  }

  implicit class StringToClobOption(val str: String) extends AnyVal {
    /*
    Many, many Clob implementations have problems with empty char arrays.
    http://hg.openjdk.java.net/jdk8/jdk8/jdk/file/687fd7c7986d/src/share/classes/javax/sql/rowset/serial/SerialClob.java#l268
    http://hg.openjdk.java.net/jdk9/client/jdk/file/1158c3e5bd9c/src/java.sql.rowset/share/classes/javax/sql/rowset/serial/SerialClob.java#l270
    https://github.com/apache/derby/blob/10.13/java/engine/org/apache/derby/iapi/types/HarmonySerialClob.java#L109
    https://github.com/arteam/hsqldb/blob/2.3.4/src/org/hsqldb/jdbc/JDBCClob.java#L196
     */

    import eu.timepit.refined.api.Refined
    import eu.timepit.refined.collection.NonEmpty

    def toClobOption: Option[SerialClob] = if (str.isEmpty) None else Option(new SerialClob(str.toCharArray))

    def toClob(default: String Refined NonEmpty): SerialClob = {
      val nonEmpty = if (str.isEmpty) default.value else str
      new SerialClob(nonEmpty.toCharArray)
    }
  }

  implicit class BlobToBytes(val blob: Blob) extends AnyVal {
    def toBytes: Array[Byte] = {
      // See notes on empty blob issues in BytesOptionToBlob
      val length = blob.length.toInt
      // yes, it starts at 1
      if (length == 0) Array.empty else blob.getBytes(1, length)
    }
  }

  implicit class BlobOptionToBytes(val blobOption: Option[Blob]) extends AnyVal {
    def toBytesOption: Option[Array[Byte]] = blobOption.map(_.toBytes)

    def toBytes: Array[Byte] = toBytesOption.getOrElse(Array.empty)
  }

  implicit class BytesOptionToBlob(val bytesOption: Option[Array[Byte]]) extends AnyVal {
    /*
    Many, many Blob implementations (but fewer than Clob) have problems with empty byte arrays.
    http://hg.openjdk.java.net/jdk8/jdk8/jdk/file/687fd7c7986d/src/share/classes/javax/sql/rowset/serial/SerialBlob.java#l178
    http://hg.openjdk.java.net/jdk9/client/jdk/file/1158c3e5bd9c/src/java.sql.rowset/share/classes/javax/sql/rowset/serial/SerialBlob.java#l178
    https://github.com/apache/derby/blob/10.13/java/engine/org/apache/derby/iapi/types/HarmonySerialBlob.java#L111
    OK! -> https://github.com/arteam/hsqldb/blob/2.3.4/src/org/hsqldb/jdbc/JDBCBlob.java#L184
     */
    def toBlobOption: Option[SerialBlob] = bytesOption.flatMap(_.toBlobOption)
  }

  implicit class BytesToBlobOption(val bytes: Array[Byte]) extends AnyVal {
    def toBlobOption: Option[SerialBlob] = if (bytes.isEmpty) None else Option(new SerialBlob(bytes))
  }

  implicit class EnhancedFiniteDuration(val duration: FiniteDuration) extends AnyVal {
    def ago: Timestamp = new Timestamp(System.currentTimeMillis() - duration.toMillis)
  }
}
