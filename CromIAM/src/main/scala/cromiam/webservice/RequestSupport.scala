package cromiam.webservice

import akka.http.scaladsl.model.HttpRequest
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import cromiam.auth.User
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId
import akka.http.scaladsl.model.HttpResponse
import akka.http.scaladsl.server.Directives.{authorize, complete, onComplete}
import akka.http.scaladsl.server.Route
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient

import scala.util.{Failure, Success}

trait RequestSupport {
  def extractStrictRequest: Directive1[HttpRequest] = {
    toStrictEntity(Timeout) tflatMap { _ =>
      extractRequest flatMap { request =>
        provide(request)
      }
    }
  }

  /**
    * Obtain both the user id header from the proxy as well as the bearer token and pass that back
    * into the route logic as a User object
    */
  def extractUser: Directive1[User] = {
    (headerValueByName("OIDC_CLAIM_user_id") & headerValuePF { case a: Authorization => a }) tmap { case (userId, auth) =>
      User(WorkbenchUserId(userId), auth)
    }
  }

  def extractUserAndStrictRequest: Directive[(User, HttpRequest)] = {
    for {
      user <- extractUser
      request <- extractStrictRequest
    } yield (user, request)
  }

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
        val message = s"Unable to look up enablement status for user ${user.userId}: ${e.getMessage}. Please try again later."
        throw new RuntimeException(message, e)
    }
  }
}
