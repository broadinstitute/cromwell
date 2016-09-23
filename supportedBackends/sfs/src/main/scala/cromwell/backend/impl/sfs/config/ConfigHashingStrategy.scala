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

  // TODO clean up this mess
  def apply(config: Config): ConfigHashingStrategy = {
    config.getConfigOption("filesystems.local.hashing") map { hashingConfig =>
      val checkSiblingMd5 = hashingConfig.getBooleanOr("check-sibling-md5", default = false)
      val checkSiblingMd5Message = if (checkSiblingMd5) {
        "Check for sibling md5 first."
      } else {
        "Do not check for sibling md5."
      }

      hashingConfig.getStringOr("strategy", "path") match {
        case "path" =>
          logger.info(s"Using path hashing strategy. $checkSiblingMd5Message")
          HashPathStrategy(checkSiblingMd5)
        case "file" =>
          logger.info(s"Using file hashing strategy. $checkSiblingMd5Message")
          HashFileStrategy(checkSiblingMd5)
        case what =>
          logger.warn(s"Unrecognized hashing strategy $what. Defaulting to path hashing.")
          HashPathStrategy(checkSiblingMd5)
      }
    } getOrElse {
      logger.info(s"Using path hashing strategy. Do not check for sibling md5.")
      HashPathStrategy(false)
    }
  }
}

sealed trait ConfigHashingStrategy {
  def checkSiblingMd5: Boolean
  def getHash(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] = {
    val file = File(request.file.valueString)

    if (checkSiblingMd5) {
      precomputedMd5(file) match {
        case Some(md5) => Try(md5.contentAsString)
        case None => hash(file)
      }
    } else hash(file)
  }

  protected def hash(file: File): Try[String]

  private def precomputedMd5(file: File): Option[File] = {
    file.siblings find { _.name == s"${file.name}.md5" }
  }
}

case class HashPathStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override def hash(file: File): Try[String] = {
    Try(org.apache.commons.codec.digest.DigestUtils.md5Hex(file.path.toAbsolutePath.toString))
  }
}
case class HashFileStrategy(checkSiblingMd5: Boolean) extends ConfigHashingStrategy {
  override protected def hash(file: File): Try[String] = {
    tryWithResource(() => file.newInputStream) { inputStream =>
      org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream)
    }
  }
}
