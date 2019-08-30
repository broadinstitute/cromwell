package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime
import java.time.{LocalDateTime, ZoneOffset}

import cats.effect.IO
import cloud.nio.spi.CloudNioRegularFileAttributes
import org.apache.commons.lang3.exception.ExceptionUtils


class DrsCloudNioRegularFileAttributes(drsPath: String, drsPathResolver: DrsPathResolver) extends CloudNioRegularFileAttributes{

  private def throwRuntimeException(missingKey: String) = {
    new RuntimeException(s"Failed to resolve DRS path $drsPath. The response from Martha doesn't contain the key '$missingKey'.")
  }

  private def convertToFileTime(timeInString: String): IO[FileTime] = {
    //Here timeInString is assumed to be a ISO-8601 DateTime without timezone
    IO(LocalDateTime.parse(timeInString).toInstant(ZoneOffset.UTC)).map(FileTime.from).handleErrorWith {
      e => IO.raiseError(new RuntimeException(s"Error while parsing 'updated' value from Martha to FileTime for DRS path $drsPath. Reason: ${ExceptionUtils.getMessage(e)}."))
    }
  }


  override def fileHash: Option[String] = {
    drsPathResolver.resolveDrsThroughMartha(drsPath).map(marthaResponse => {
      marthaResponse.dos.data_object.checksums.flatMap {
        _.collectFirst{ case c if c.`type`.equalsIgnoreCase("md5") => c.checksum }
      }
    }).unsafeRunSync()
  }


  override def lastModifiedTime(): FileTime = {
    val lastModifiedIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      lastModifiedInString <- IO.fromEither(marthaResponse.dos.data_object.updated.toRight(throwRuntimeException("updated")))
      lastModified <- convertToFileTime(lastModifiedInString)
    } yield lastModified

    lastModifiedIO.unsafeRunSync()
  }


  override def size(): Long = {
    val sizeIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      size <- IO.fromEither(marthaResponse.dos.data_object.size.toRight(throwRuntimeException("size")))
    } yield size

    sizeIO.unsafeRunSync()
  }


  override def fileKey(): String = drsPath
}
