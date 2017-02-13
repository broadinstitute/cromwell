package cromwell.backend.impl.tes

import akka.event.LoggingAdapter
import cromwell.backend.callcaching.FileHashingActor.SingleFileHashRequest
import cromwell.core.path.DefaultPathBuilder
import cromwell.util.TryWithResource._

import scala.language.postfixOps
import scala.util.Try

private[tes] object TesBackendFileHashing {
  def getMd5Result(request: SingleFileHashRequest, log: LoggingAdapter): Try[String] = {
    val path = DefaultPathBuilder.build(request.file.valueString) recover {
      case failure => throw new RuntimeException("Failed to construct path to hash", failure)
    } get

    tryWithResource(() => path.newInputStream) { inputStream =>
      org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream)
    }
  }
}
