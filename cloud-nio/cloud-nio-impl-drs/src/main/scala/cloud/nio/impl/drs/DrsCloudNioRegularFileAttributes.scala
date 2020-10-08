package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime
import java.time.{LocalDateTime, OffsetDateTime, ZoneOffset}

import cats.effect.IO
import cloud.nio.spi.CloudNioRegularFileAttributes
import org.apache.commons.lang3.exception.ExceptionUtils

class DrsCloudNioRegularFileAttributes(drsPath: String,
                                       sizeOption: Option[Long],
                                       hashOption: Option[String],
                                       timeCreatedOption: Option[FileTime],
                                       timeUpdatedOption: Option[FileTime],
                                      ) extends CloudNioRegularFileAttributes{

  override def fileKey(): String = drsPath

  override def size(): Long = sizeOption.getOrElse(0)

  override def fileHash: Option[String] = hashOption

  override def creationTime(): FileTime = timeCreatedOption.getOrElse(lastModifiedTime())

  override def lastModifiedTime(): FileTime = timeUpdatedOption.getOrElse(FileTime.fromMillis(0))
}

object DrsCloudNioRegularFileAttributes {
  private val priorityHashList: Seq[String] = Seq("crc32c", "md5", "sha256")

  def getPreferredHash(hashesOption: Option[Map[String, String]]): Option[String] = {
    hashesOption match {
      case Some(hashes) if hashes.nonEmpty =>
        val preferredHash = priorityHashList collectFirst {
          case hashKey if hashes.contains(hashKey) => hashes(hashKey)
        }

        // if no preferred hash was found, sort the hashes alphabetically by type and take the first one
        Option(
          preferredHash.getOrElse(
            hashes.toSeq minBy {
              case (hashType, _) => hashType
            } match {
              case (_, hashValue) => hashValue
            }
          )
        )
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
