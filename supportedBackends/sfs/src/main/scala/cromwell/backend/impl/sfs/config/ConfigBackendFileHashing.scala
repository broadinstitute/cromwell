package cromwell.backend.impl.sfs.config

import java.io.{File, FileInputStream}

import akka.event.LoggingAdapter
import cromwell.backend.callcaching.FileHasherWorkerActor.SingleFileHashRequest

import scala.util.Try
import cromwell.util.TryWithResource._

private[config] object ConfigBackendFileHashing {
  def getMd5Result(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] =
    tryWithResource(() => new FileInputStream(new File(request.file.valueString))) { fileInputStream =>
      org.apache.commons.codec.digest.DigestUtils.md5Hex(fileInputStream)
    }
}
