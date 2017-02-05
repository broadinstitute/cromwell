package cromwell.backend.impl.jes

import com.google.api.client.http.HttpResponseException
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.core.path.Path

package object io {
  implicit class PathEnhanced(val path: Path) extends AnyVal {

    def writeAsJson(content: String): Path = {
      path.writeBytes(content.getBytes.toIterator)(Seq(CloudStorageOptions.withMimeType("application/json")))
    }

    def writeAsText(content: String): Path = {
      path.writeBytes(content.getBytes.toIterator)(Seq(CloudStorageOptions.withMimeType("text/plain")))
    }
  }

  private [jes] def isFatalJesException(t: Throwable): Boolean = t match {
    case e: HttpResponseException if e.getStatusCode == 403 => true
    case e: HttpResponseException if e.getStatusCode == 400 && e.getContent.contains("INVALID_ARGUMENT") => true
    case _ => false
  }

  private [jes] def isTransientJesException(t: Throwable): Boolean = t match {
    // Quota exceeded
    case e: HttpResponseException if e.getStatusCode == 429 => true
    case _ => false
  }
}
