package cromwell.backend.impl.sfs.config

import akka.event.LoggingAdapter
import better.files._
import cromwell.backend.callcaching.FileHashingActor.SingleFileHashRequest
import cromwell.util.TryWithResource._

import scala.util.Try

private[config] object ConfigBackendFileHashing {
  def getMd5Result(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] =
    tryWithResource(() => File(request.file.valueString).newInputStream) { inputStream =>
      org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream)
    }
}
