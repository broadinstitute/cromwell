package cromwell.backend.impl.sfs.config

import akka.event.LoggingAdapter
import better.files._
import com.typesafe.config.Config
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.callcaching.FileHashingActor.SingleFileHashRequest
import cromwell.util.TryWithResource._

import scala.util.Try

private[config] object ConfigBackendFileHashing {
  def getMd5Result(config: Config)(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] = {
    val file = File(request.file.valueString)
    precomputedMd5(file) match {
      case Some(md5) => md5.contentAsString
      case None => config.g
    }
    tryWithResource(() => File(request.file.valueString).newInputStream) { inputStream =>
      org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream)
    }
  }

  private def precomputedMd5(file: File): Option[File] = {
    file.siblings find { _.name == s"${file.name}.md5" }
  }
}
