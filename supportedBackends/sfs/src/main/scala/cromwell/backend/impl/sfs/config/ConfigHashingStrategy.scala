package cromwell.backend.impl.sfs.config

import akka.event.LoggingAdapter
import better.files.File
import com.typesafe.config.Config
import cromwell.backend.callcaching.FileHashingActor.SingleFileHashRequest
import cromwell.util.TryWithResource._
import lenthall.config.ScalaConfig._
import org.slf4j.LoggerFactory

import scala.util.Try

object ConfigHashingStrategy {
  val logger = LoggerFactory.getLogger(getClass)
  val defaultStrategy = HashPathStrategy(false)

  def apply(hashingConfig: Config): ConfigHashingStrategy = {
      val checkSiblingMd5 = hashingConfig.getBooleanOr("check-sibling-md5", default = false)

      hashingConfig.getStringOr("strategy", "path") match {
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
    val file = File(request.file.valueString)

    if (checkSiblingMd5) {
      precomputedMd5(file) match {
        case Some(md5) => Try(md5.contentAsString)
        case None => hash(file)
      }
    } else hash(file)
  }


  private def precomputedMd5(file: File): Option[File] = {
    file.siblings find { _.name == s"${file.name}.md5" }
  }

  override def toString = {
    s"Call caching hashing strategy: $checkSiblingMessage$description."
  }
}

case class HashPathStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override def hash(file: File): Try[String] = {
    Try(org.apache.commons.codec.digest.DigestUtils.md5Hex(file.path.toAbsolutePath.toString))
  }

  override val description = "hash file path"
}
case class HashFileStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override protected def hash(file: File): Try[String] = {
    tryWithResource(() => file.newInputStream) { inputStream =>
      org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream)
    }
  }

  override val description = "hash file content"
}
