package cromiam.webservice

import akka.http.scaladsl.model.HttpRequest
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import cromiam.auth.User
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId

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

  def extractUserAndRequest: Directive[(User, HttpRequest)] = {
    for {
      user <- extractUser
      request <- extractStrictRequest
    } yield (user, request)
  }
}
