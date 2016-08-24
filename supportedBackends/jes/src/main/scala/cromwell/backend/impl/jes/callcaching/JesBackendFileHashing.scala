package cromwell.backend.impl.jes.callcaching

import akka.event.LoggingAdapter
import cromwell.backend.callcaching.FileHashingActor.SingleFileHashRequest
import cromwell.backend.impl.jes.JesBackendInitializationData

import scala.util.{Failure, Try}

private[jes] object JesBackendFileHashing {
  def getCrc32c(singleFileHashRequest: SingleFileHashRequest, log: LoggingAdapter): Try[String] = {
    def usingJesInitData(jesInitData: JesBackendInitializationData) = for {
      path <- Try(jesInitData.workflowPaths.gcsFileSystemWithUserAuth.getPath(singleFileHashRequest.file.valueString))
      crc32c <- Try(jesInitData.workflowPaths.gcsFileSystemWithUserAuth.gcsFileSystemProvider.crc32cHash(path))
    } yield crc32c

    singleFileHashRequest.initializationData match {
      case Some(jesInitData: JesBackendInitializationData) => usingJesInitData(jesInitData)
      case _ => Failure(new IllegalArgumentException("Need JesBackendInitializationData to generate a GCS CRC32C hash"))
    }
  }
}
