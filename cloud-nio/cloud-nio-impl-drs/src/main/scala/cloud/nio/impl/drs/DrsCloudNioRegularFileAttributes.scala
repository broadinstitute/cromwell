package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime
import java.time.{LocalDateTime, OffsetDateTime, ZoneOffset}
import cats.effect.IO
import cloud.nio.spi.{CloudNioRegularFileAttributes, FileHash}
import org.apache.commons.lang3.exception.ExceptionUtils

class DrsCloudNioRegularFileAttributes(drsPath: String,
                                       sizeOption: Option[Long],
                                       hashOption: Option[FileHash],
                                       timeCreatedOption: Option[FileTime],
                                       timeUpdatedOption: Option[FileTime],
                                      ) extends CloudNioRegularFileAttributes{

  override def fileKey(): String = drsPath

  override def size(): Long = sizeOption.getOrElse(0)

  override def fileHash: Option[FileHash] = hashOption

//  override def hashType: Option[String] = hashTypeOption

  override def creationTime(): FileTime = timeCreatedOption.getOrElse(lastModifiedTime())

  override def lastModifiedTime(): FileTime = timeUpdatedOption.getOrElse(FileTime.fromMillis(0))
}

object DrsCloudNioRegularFileAttributes {
  private val priorityHashList: Seq[FileHash.Value] = Seq(FileHash.Crc32c, FileHash.Md5, FileHash.Sha256, FileHash.Etag)

  def getPreferredHash(hashesOption: Option[Map[String, String]]): Option[FileHash] = {
    hashesOption match {
      case Some(hashes: Map[String, String]) if hashes.nonEmpty =>
        priorityHashList collectFirst {
          case hashKey if hashes.contains(hashKey.toString) => FileHash(hashKey, hashes(hashKey.toString))
        }
        // if no preferred hash was found, go ahead and return none because we don't support anything that the DRS object is offering
      case _ => None
    }
  }

  private def convertToOffsetDateTime(timeInString: String): IO[OffsetDateTime] = {
    // Here timeInString is assumed to be a ISO-8601 DateTime with timezone
    IO(OffsetDateTime.parse(timeInString))
      .handleErrorWith(
        offsetDateTimeException =>
          // As a fallback timeInString is assumed to be a ISO-8601 DateTime without timezone
          IO(LocalDateTime.parse(timeInString).atOffset(ZoneOffset.UTC))
            .handleErrorWith(_ => IO.raiseError(offsetDateTimeException))
      )
  }

  private def convertToFileTime(timeInString: String): IO[FileTime] = {
    convertToOffsetDateTime(timeInString)
      .map(_.toInstant)
      .map(FileTime.from)
  }

  def convertToFileTime(drsPath: String, key: MarthaField.Value, timeInStringOption: Option[String]): IO[Option[FileTime]] = {
    timeInStringOption match {
      case None => IO.pure(None)
      case Some(timeInString) =>
        convertToFileTime(timeInString)
          .map(Option(_))
          .handleErrorWith(
            throwable =>
              IO.raiseError(
                new RuntimeException(
                  s"Error while parsing '$key' value from Martha to FileTime for DRS path $drsPath. " +
                    s"Reason: ${ExceptionUtils.getMessage(throwable)}.",
                  throwable,
                )
              )
          )
    }
  }
}
