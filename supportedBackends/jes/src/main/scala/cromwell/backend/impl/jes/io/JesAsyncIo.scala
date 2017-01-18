package cromwell.backend.impl.jes.io

import java.nio.file.Path

import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.services.io.AsyncIo

import scala.concurrent.Future

trait JesAsyncIo extends AsyncIo {
  def writeAsJson(file: Path, content: String): Future[Unit] = {
    write(file, content, Seq(CloudStorageOptions.withMimeType("application/json")))
  }
}
