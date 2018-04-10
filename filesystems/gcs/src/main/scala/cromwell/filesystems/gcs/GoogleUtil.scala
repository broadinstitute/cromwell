package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.client.http.HttpResponseException
import com.google.auth.Credentials
import cromwell.cloudsupport.gcp.auth.{GoogleAuthMode, OptionLookupException}
import cromwell.core.retry.Retry
import cromwell.core.{CromwellFatalExceptionMarker, WorkflowOptions}

import scala.concurrent.{ExecutionContext, Future}

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

  implicit class EnhancedGoogleAuthMode(val googleAuthMode: GoogleAuthMode) extends AnyVal {
    /**
      * Retries getting the credentials three times.
      */
    def retryCredential(options: WorkflowOptions)
                       (implicit as: ActorSystem, ec: ExecutionContext): Future[Credentials] = {
      def credential(): Credentials = {
        try {
          googleAuthMode.credential((key: String) => options.get(key).get)
        } catch {
          case exception: OptionLookupException =>
            throw new IllegalArgumentException(s"Missing parameters in workflow options: ${exception.key}", exception)
              with CromwellFatalExceptionMarker
        }
      }

      def isFatal(throwable: Throwable): Boolean = {
        throwable match {
          case _: IllegalArgumentException => Option(throwable.getCause).exists(isFatal)
          case _ => GoogleAuthMode.isFatal(throwable)
        }
      }

      Retry.withRetry(
        () => Future(credential()),
        isFatal = isFatal,
        maxRetries = Option(3)
      )
    }
  }
}
