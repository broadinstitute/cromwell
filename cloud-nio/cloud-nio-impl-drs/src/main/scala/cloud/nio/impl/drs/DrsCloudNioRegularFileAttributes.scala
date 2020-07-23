package cloud.nio.impl.drs

import java.nio.file.attribute.FileTime
import java.time.{LocalDateTime, ZoneOffset}

import cats.effect.IO
import cloud.nio.spi.CloudNioRegularFileAttributes
import org.apache.commons.lang3.exception.ExceptionUtils


class DrsCloudNioRegularFileAttributes(drsPath: String, drsPathResolver: DrsPathResolver) extends CloudNioRegularFileAttributes{

  private def createMissingKeyException(missingKey: String) = {
    new RuntimeException(s"Failed to resolve DRS path $drsPath. The response from Martha doesn't contain the key '$missingKey'.")
  }

  private def convertToFileTime(timeInString: String): IO[FileTime] = {
    //Here timeInString is assumed to be a ISO-8601 DateTime without timezone
    IO(LocalDateTime.parse(timeInString).toInstant(ZoneOffset.UTC)).map(FileTime.from).handleErrorWith {
      e => IO.raiseError(new RuntimeException(s"Error while parsing 'updated' value from Martha to FileTime for DRS path $drsPath. Reason: ${ExceptionUtils.getMessage(e)}."))
    }
  }


  override def fileHash: Option[String] = {
    /* ###### ANOTHER APPROACH FOR GETTING HASHES #######
    val predicates = List[String => Boolean](
      x => x.equalsIgnoreCase("crc32c"),
      x => x.equalsIgnoreCase("md5"),
      x => x.equalsIgnoreCase("sha256")
    )
    val indexedPredicates = predicates.reverse.zipWithIndex
    def score(x: String): Option[Int] = indexedPredicates.find(_._1(x)).map(_._2)

    def priorityFind( l: List[String] ): Option[String] = {
      val filtered = l.view.flatMap{x => score(x).map(x -> _) }
      if ( filtered.isEmpty ) None
      else Some( filtered.maxBy(_._2)._1 )
    }
     */

    drsPathResolver.resolveDrsThroughMartha(drsPath).map((marthaResponse: MarthaResponse) => {

      println(s"############## MARTHA RESPONSE ##############")
      println(marthaResponse.hashes)

      val maybeHashes: Option[Map[String, String]] = marthaResponse.hashes
      val hashes: Map[String, String] = maybeHashes match {
        case Some(value) => value
        case None => throw createMissingKeyException("hashes")
      }

      val preferredHash: Option[String] = List("crc32c", "md5", "sha256").collectFirst {
        case k if hashes.contains(k) => hashes(k)
      }

      preferredHash match {
        case Some(hash) =>
          Option(hash)
        case None =>
          Option(hashes.toSeq.minBy(_._1)._2)
      }
    }).unsafeRunSync()
  }


  override def lastModifiedTime(): FileTime = {
    val lastModifiedIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      lastModifiedInString <- IO.fromEither(marthaResponse.timeUpdated.toRight(createMissingKeyException("updated")))
      lastModified <- convertToFileTime(lastModifiedInString)
    } yield lastModified

    lastModifiedIO.unsafeRunSync()
  }


  override def size(): Long = {
    val sizeIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      size <- IO.fromEither(marthaResponse.size.toRight(createMissingKeyException("size")))
    } yield size

    sizeIO.unsafeRunSync()
  }


  override def fileKey(): String = drsPath
}
