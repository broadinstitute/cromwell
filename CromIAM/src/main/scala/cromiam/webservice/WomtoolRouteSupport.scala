package cromiam.webservice

import akka.http.scaladsl.server.Directives._
import cromiam.cromwell.CromwellClient
import cromwell.api.model._

trait WomtoolRouteSupport extends RequestSupport {
  val cromwellClient: CromwellClient

  val womtoolRoutes =
    path("api" / "womtool" / Segment / "describe") { _ =>
      post {
        extractStrictRequest { req =>
          complete {
            // My understanding is that the proxies enforce authentication for any /api route, and that no further
            // authentication nor authorization is required here since any authenticated user can describe workflows.
            cromwellClient.forwardToCromwell(req).asHttpResponse
          }
        }
      }
    }

}
