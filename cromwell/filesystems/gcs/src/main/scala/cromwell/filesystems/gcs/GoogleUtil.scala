package cromwell.filesystems.gcs

import com.google.api.client.http.HttpResponseException

object GoogleUtil {
  /**
    * Extract status code from an exception if it's a com.google.api.client.http.HttpResponseException
    */
  def extractStatusCode(exception: Throwable): Option[Int] = {
    exception match {
      case t: HttpResponseException => Option(t.getStatusCode)
      case _ => None
    }
  }
}
