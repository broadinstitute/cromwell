package cromwell.backend.google.batch

import com.google.api.client.http.HttpResponseException
import com.google.api.gax.rpc.ApiException
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.core.path.Path

package object io {
  implicit class PathEnhanced(val path: Path) extends AnyVal {

    def writeAsJson(content: String): Path =
      path.writeBytes(content.getBytes.iterator)(Seq(CloudStorageOptions.withMimeType("application/json")))

    def writeAsText(content: String): Path =
      path.writeBytes(content.getBytes.iterator)(Seq(CloudStorageOptions.withMimeType("text/plain")))
  }

  private[batch] def isFatalBatchException(t: Throwable): Boolean = t match {
    case e: ApiException if e.getStatusCode.getCode.getHttpStatusCode == 403 => true
    case e: HttpResponseException if e.getStatusCode == 403 => true
    case e: HttpResponseException if e.getStatusCode == 400 && e.getContent.contains("INVALID_ARGUMENT") => true
    case _ => false
  }

  private[batch] def isTransientBatchException(t: Throwable): Boolean = t match {
    // Quota exceeded
    case e: HttpResponseException if e.getStatusCode == 429 => true
    case e: ApiException if e.getStatusCode.getCode.getHttpStatusCode == 429 => true
    case _ => false
  }
}
