package cromwell.backend.impl.sfs.config

import akka.event.LoggingAdapter
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.core.path.DefaultPathBuilder
import cromwell.util.TryWithResource._

import scala.language.postfixOps
import scala.util.Try

private[config] object ConfigBackendFileHashing {
  def getMd5Result(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] = {
    val path = DefaultPathBuilder.build(request.file.valueString) recover { case failure =>
      throw new RuntimeException("Failed to construct path to hash", failure)
    } get

    tryWithResource(() => path.newInputStream) { inputStream =>
      org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream)
    }
  }

}
