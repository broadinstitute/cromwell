package cromwell.filesystems.gcs

import akka.actor.ActorSystem
import com.google.api.client.http.HttpResponseException
import com.google.auth.Credentials
import com.google.cloud.BaseServiceException
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
      case t: BaseServiceException => Option(t.getCode)
      case _ => None
    }
  }

  implicit class EnhancedGoogleAuthMode(val googleAuthMode: GoogleAuthMode) extends AnyVal {
    /**
      * Retries getting the credentials three times.
      *
      * There is nothing GCS specific about this method. This package just happens to be the lowest level with access
      * to core's version of Retry + cloudSupport's implementation of GoogleAuthMode.
      */
    def retryCredentials(options: WorkflowOptions, scopes: Iterable[String])
                        (implicit actorSystem: ActorSystem, executionContext: ExecutionContext): Future[Credentials] = {
      def credential(): Credentials = {
        try {
          googleAuthMode.credentials(options.get(_).get, scopes)
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
