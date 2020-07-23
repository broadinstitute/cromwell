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

    drsPathResolver.resolveDrsThroughMartha(drsPath).map(marthaResponse => {

      println(s"############## MARTHA RESPONSE ##############")
      println(marthaResponse.hashes)

      marthaResponse.hashes.flatMap { h =>

        val hashMap: Map[String, Json] = h.toMap

//        val prefHash = priorityFind(hashMap.keys.toList)
//        println(prefHash)

        val preferredHash: Option[Json] = hashMap.keys.toList match {
          case k if k.contains("crc32c") => hashMap.get("crc32c")
          case k if k.contains("md5") => hashMap.get("md5")
          case k if k.contains("sha256") => hashMap.get("sha256")
          case _ => None
        }

//        val preferredHash: Option[Json] = hashMap.keys.collectFirst {
//          case k if k.equalsIgnoreCase("md5") =>
//            println("hit md5")
//            hashMap.get(k)
//          case k if k.equalsIgnoreCase("crc32c") =>
//            println("hit crc")
//            hashMap.get(k)
//          case k if k.equalsIgnoreCase("sha256") =>
//            println("hit sha")
//            hashMap.get(k)
//        }.flatten

        preferredHash match {
          case Some(hash) =>
            hash.asString
          case None =>
            hashMap.toSeq.minBy(_._1)._2.asString
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
