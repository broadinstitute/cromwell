package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime
import java.time.{LocalDateTime, ZoneOffset}

import cats.effect.IO
import cloud.nio.spi.CloudNioRegularFileAttributes
import com.typesafe.config.Config


class DrsCloudNioRegularFileAttributes(config: Config, drsPath: String, drsPathResolver: DrsPathResolver) extends CloudNioRegularFileAttributes{

  private def throwRuntimeException(missingKey: String) = {
    new RuntimeException(s"Failed to resolve DRS path $drsPath. The response from Martha doesn't contain the key '$missingKey'.")
  }

  private def timeInStringToFileTime(timeInString: String): FileTime = {
    //Here timeInString is assumed to be a ISO-8601 DateTime without timezone
    val instant = LocalDateTime.parse(timeInString).toInstant(ZoneOffset.UTC)
    FileTime.from(instant)
  }


  override def fileHash: Option[String] = {
    drsPathResolver.resolveDrsThroughMartha(drsPath).map(marthaResponse => {
      val checksumsArrayOption = marthaResponse.dos.data_object.checksums

      checksumsArrayOption match {
        case Some(checksumsArray) => checksumsArray.collectFirst{ case c if c.`type`.equalsIgnoreCase("md5") => c.checksum }
        case None => None
      }
    }).unsafeRunSync()
  }


  override def lastModifiedTime(): FileTime = {
    val lastModifiedInStringIO =  for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      lastModifiedInString <- IO.fromEither(marthaResponse.dos.data_object.updated.map(timeInStringToFileTime).toRight(throwRuntimeException("updated")))
    } yield lastModifiedInString

    lastModifiedInStringIO.unsafeRunSync()
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
