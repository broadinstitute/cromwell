package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime
import java.time.{LocalDateTime, ZoneOffset}

import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioRegularFileAttributes._
import cloud.nio.spi.CloudNioRegularFileAttributes
import org.apache.commons.lang3.exception.ExceptionUtils


class DrsCloudNioRegularFileAttributes(drsPath: String, drsPathResolver: EngineDrsPathResolver) extends CloudNioRegularFileAttributes{

  private def convertToFileTime(timeInString: String): IO[FileTime] = {
    //Here timeInString is assumed to be a ISO-8601 DateTime without timezone
    IO(LocalDateTime.parse(timeInString).toInstant(ZoneOffset.UTC)).map(FileTime.from).handleErrorWith {
      e => IO.raiseError(new RuntimeException(s"Error while parsing 'updated' value from Martha to FileTime for DRS path $drsPath. Reason: ${ExceptionUtils.getMessage(e)}."))
    }
  }


  override def fileHash: Option[String] = {
    drsPathResolver.resolveDrsThroughMartha(drsPath).map( marthaResponse => {
      marthaResponse.hashes match {
        case Some(hashes) => getPreferredHash(hashes)
        case None => throw createMissingKeyException(drsPath, "hashes")
      }
    }).unsafeRunSync()
  }


  override def lastModifiedTime(): FileTime = {
    val lastModifiedIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      lastModifiedInString <- IO.fromEither(marthaResponse.timeUpdated.toRight(createMissingKeyException(drsPath, "updated")))
      lastModified <- convertToFileTime(lastModifiedInString)
    } yield lastModified

    lastModifiedIO.unsafeRunSync()
  }


  override def size(): Long = {
    val sizeIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      size <- IO.fromEither(marthaResponse.size.toRight(createMissingKeyException(drsPath, "size")))
    } yield size

    sizeIO.unsafeRunSync()
  }


  override def fileKey(): String = drsPath
}

object DrsCloudNioRegularFileAttributes {
  private val priorityHashList: Seq[String] = Seq("crc32c", "md5", "sha256")

  def getPreferredHash(hashes: Map[String, String]): Option[String] = {
    val preferredHash: Option[String] = priorityHashList.collectFirst {
      case hashKey if hashes.contains(hashKey) => hashes(hashKey)
    }

    // if no preferred hash was found, sort the hashes alphabetically by type and take the first one
    Option(preferredHash.getOrElse(hashes.toSeq.minBy(_._1)._2))
  }

  def createMissingKeyException(drsPath: String, missingKey: String) = {
    new RuntimeException(s"Failed to resolve DRS path $drsPath. The response from Martha doesn't contain the key '$missingKey'.")
  }
}

