package cromwell.backend.impl.sfs.config

import akka.event.LoggingAdapter
import better.files.File
import com.typesafe.config.Config
import cromwell.backend.callcaching.FileHashingActor.SingleFileHashRequest
import cromwell.backend.sfs.SharedFileSystemBackendInitializationData
import cromwell.core.path.PathFactory
import cromwell.util.TryWithResource._
import cromwell.util.FileUtil._
import net.ceedubs.ficus.Ficus._
import org.apache.commons.codec.digest.DigestUtils
import org.slf4j.LoggerFactory

import scala.util.{Failure, Try}

object ConfigHashingStrategy {
  val logger = LoggerFactory.getLogger(getClass)
  val defaultStrategy = HashFileStrategy(false)

  def apply(hashingConfig: Config): ConfigHashingStrategy = {
      val checkSiblingMd5 = hashingConfig.as[Option[Boolean]]("check-sibling-md5").getOrElse(false)

      hashingConfig.as[Option[String]]("hashing-strategy").getOrElse("file") match {
        case "path" => HashPathStrategy(checkSiblingMd5)
        case "file" => HashFileStrategy(checkSiblingMd5)
        case what =>
          logger.warn(s"Unrecognized hashing strategy $what.")
          HashPathStrategy(checkSiblingMd5)
      }
  }
}

abstract class ConfigHashingStrategy {
  def checkSiblingMd5: Boolean
  protected def hash(file: File): Try[String]
  protected def description: String

  protected lazy val checkSiblingMessage = if (checkSiblingMd5) "Check first for sibling md5 and if not found " else ""

  def getHash(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] = {
    def usingSFSInitData(initData: SharedFileSystemBackendInitializationData) = {
      val pathBuilders = initData.workflowPaths.pathBuilders
      val file = PathFactory.buildFile(request.file.valueString, pathBuilders).followSymlinks

      if (checkSiblingMd5) {
        precomputedMd5(file) match {
          case Some(md5) => Try(md5.contentAsString)
          case None => hash(file)
        }
      } else hash(file)
    }

    request.initializationData match {
      case Some(initData: SharedFileSystemBackendInitializationData) => usingSFSInitData(initData)
      case _ => Failure(new IllegalArgumentException("Need SharedFileSystemBackendInitializationData to calculate hash."))
    }
  }

  private def precomputedMd5(file: File): Option[File] = {
    val md5 = file.sibling(s"${file.name}.md5")
    if (md5.exists) Option(md5) else None
  }

  override def toString = {
    s"Call caching hashing strategy: $checkSiblingMessage$description."
  }
}

final case class HashPathStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override def hash(file: File): Try[String] = {
    Try(DigestUtils.md5Hex(file.path.toAbsolutePath.toString))
  }

  override val description = "hash file path"
}

final case class HashFileStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override protected def hash(file: File): Try[String] = {
    tryWithResource(() => file.newInputStream) { DigestUtils.md5Hex }
  }

  override val description = "hash file content"
}
