package cromiam.webservice

import akka.http.scaladsl.model.{HttpResponse, StatusCodes}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient.SamDenialException
import cromiam.server.config.CromIamServerConfig
import cromwell.api.model._

trait AdminRouteSupport extends RequestSupport {
  // When this trait is mixed into `CromIamApiService` the value of `cromwellClient` is one of the readers
  val cromwellClient: CromwellClient

  protected def configuration: CromIamServerConfig

  def authenticateAdminActionAndForward(): Route = extractUserAndRequest { (user, req) =>

    if (configuration.adminWhitelist.contains(user.userId)) {
      complete {
        cromwellClient.forwardToCromwell(req).asHttpResponse
      }
    } else {
      println(s"Invalid access attempt from: ${user.userId}")
      complete(HttpResponse(status = StatusCodes.Forbidden, entity = new SamDenialException().getMessage))
    }
  }

  val getAdminRoutes = List("listSubmissions")
  val postAdminRoutes = List("pauseSubmission", "pauseAll", "releaseHoldAcrossSubmission", "insertTerminalStatusInMetadata")

  val adminRoutes = concat(

    path("api"/ "admin" / Segment / "listSubmissions") { _ =>
      get {
        authenticateAdminActionAndForward()
      }
    },
    path("api"/ "admin" / Segment / "pauseSubmission") { _ =>
      post {
        authenticateAdminActionAndForward()
      }
    },
    path("api"/ "admin" / Segment / "pauseAll") { _ =>
      post {
        authenticateAdminActionAndForward()
      }
    },
    path("api"/ "admin" / Segment / "releaseHoldAcrossSubmission") { _ =>
      post {
        authenticateAdminActionAndForward()
      }
    },
    path("api"/ "admin" / Segment / "insertTerminalStatusInMetadata") { _ =>
      post {
        authenticateAdminActionAndForward()
      }
    }
  )
}
