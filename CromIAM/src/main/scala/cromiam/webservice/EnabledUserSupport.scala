package cromiam.webservice

import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.http.scaladsl.server.Directives.{authorize, complete, onComplete}
import akka.http.scaladsl.server.Route
import cromiam.auth.User
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient

import scala.util.{Failure, Success}

trait EnabledUserSupport {

  def forwardIfUserEnabled(user: User, req: HttpRequest, cromwellClient: CromwellClient, samClient: SamClient): Route = {
    import cromwell.api.model.EnhancedFailureResponseOrHttpResponseT

    onComplete(samClient.isUserEnabledSam(user, req).value.unsafeToFuture()) {
      case Success(Left(httpResponse: HttpResponse)) => complete(httpResponse)
      case Success(Right(isEnabled: Boolean)) =>
        authorize(isEnabled) {
          complete {
            cromwellClient.forwardToCromwell(req).asHttpResponse
          }
        }
      case Failure(e) =>
        val message = s"Unable to look up enablement status for user ${user.userId}: ${e.getMessage}"
        throw new RuntimeException(message, e)
    }
  }

}
