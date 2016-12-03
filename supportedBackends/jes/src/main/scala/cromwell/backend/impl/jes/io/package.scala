package cromwell.backend.impl.jes

import java.nio.file.{Files, Path}

import com.google.api.client.http.HttpResponseException
import com.google.cloud.storage.contrib.nio.CloudStorageOptions

package object io {
  implicit class PathEnhanced(val path: Path) extends AnyVal {
    import better.files._

    def writeAsJson(content: String): File = {
      Files.write(path, content.getBytes, CloudStorageOptions.withMimeType("application/json"))
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
