package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime
import java.time.{LocalDateTime, ZoneOffset}

import cats.effect.IO
import cloud.nio.spi.CloudNioRegularFileAttributes
import io.circe.Json
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
      marthaResponse.hashes.flatMap { h =>

        val hashMap: Map[String, Json] = h.toMap
        val preferredHash: Option[Json] = hashMap.keys.collectFirst {
          case k if k.equalsIgnoreCase("crc32c") => hashMap.get(k)
          case k if k.equalsIgnoreCase("md5") => hashMap.get(k)
          case k if k.equalsIgnoreCase("sha256") => hashMap.get(k)
        }.flatten

        preferredHash match {
          case Some(hash) => hash.asString
          case None => hashMap.toSeq.minBy(_._1)._2.asString
        }
      }
    }).unsafeRunSync()
  }


  override def lastModifiedTime(): FileTime = {
    val lastModifiedIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      lastModifiedInString <- IO.fromEither(marthaResponse.timeUpdated.toRight(throwRuntimeException("updated")))
      lastModified <- convertToFileTime(lastModifiedInString)
    } yield lastModified

    lastModifiedIO.unsafeRunSync()
  }


  override def size(): Long = {
    val sizeIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      size <- IO.fromEither(marthaResponse.size.toRight(throwRuntimeException("size")))
    } yield size

    sizeIO.unsafeRunSync()
  }


  override def fileKey(): String = drsPath
}
