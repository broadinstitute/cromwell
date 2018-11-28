package cromwell.backend.impl.sfs.config

import java.io.FileNotFoundException

import akka.event.LoggingAdapter
import com.typesafe.config.Config
import cromwell.backend.standard.StandardInitializationData
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.core.path.{Path, PathFactory}
import cromwell.util.TryWithResource._
import net.ceedubs.ficus.Ficus._
import org.apache.commons.codec.digest.DigestUtils
import org.slf4j.{Logger, LoggerFactory}

import scala.util.{Failure, Try}

object ConfigHashingStrategy {
  val logger: Logger = LoggerFactory.getLogger(getClass)
  val defaultStrategy = HashFileStrategy(checkSiblingMd5 = false)

  def apply(hashingConfig: Config): ConfigHashingStrategy = {
      val checkSiblingMd5 = hashingConfig.as[Option[Boolean]]("check-sibling-md5").getOrElse(false)

      hashingConfig.as[Option[String]]("hashing-strategy").getOrElse("file") match {
        case "path" => HashPathStrategy(checkSiblingMd5)
        case "file" => HashFileStrategy(checkSiblingMd5)
        case "path+modtime" => HashPathModTimeStrategy(checkSiblingMd5)
        case what =>
          logger.warn(s"Unrecognized hashing strategy $what.")
          HashPathStrategy(checkSiblingMd5)
      }
  }
}

abstract class ConfigHashingStrategy {
  def checkSiblingMd5: Boolean
  protected def hash(file: Path): Try[String]
  protected def description: String

  protected lazy val checkSiblingMessage: String =
    if (checkSiblingMd5) "Check first for sibling md5 and if not found " else ""

  def getHash(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] = {
    def usingStandardInitData(initData: StandardInitializationData) = {
      val pathBuilders = initData.workflowPaths.pathBuilders
      val file = PathFactory.buildPath(request.file.valueString, pathBuilders).followSymbolicLinks
      if (!file.exists) Failure(new FileNotFoundException(s"Cannot hash file $file because it can't be found")) else {
        if (checkSiblingMd5) {
          precomputedMd5(file) match {
            case Some(md5) => Try(md5.contentAsString.trim)
            case None => hash(file)
          }
        } else hash(file)
      }
    }

    request.initializationData match {
      case Some(initData: StandardInitializationData) => usingStandardInitData(initData)
      case _ => Failure(new IllegalArgumentException("Need SharedFileSystemBackendInitializationData to calculate hash."))
    }
  }

  private def precomputedMd5(file: Path): Option[Path] = {
    val md5 = file.sibling(s"${file.name}.md5")
    if (md5.exists) Option(md5) else None
  }

  override def toString: String = {
    s"Call caching hashing strategy: $checkSiblingMessage$description."
  }
}

final case class HashPathStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override def hash(file: Path): Try[String] = {
    Try(DigestUtils.md5Hex(file.toAbsolutePath.pathAsString))
  }

  override val description = "hash file path"
}

final case class HashPathModTimeStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override def hash(file: Path): Try[String] = {
    // Add the last modified date here to make sure these are the files we are looking for.
    Try(DigestUtils.md5Hex(file.toAbsolutePath.pathAsString + file.lastModifiedTime.toString))
  }

  override val description = "hash file path and last modified time"
}

final case class HashFileStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override protected def hash(file: Path): Try[String] = {
    tryWithResource(() => file.newInputStream) { DigestUtils.md5Hex }
  }

  override val description = "hash file content"
}
