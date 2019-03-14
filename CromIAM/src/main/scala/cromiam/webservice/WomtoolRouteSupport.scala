package cromiam.webservice

import akka.http.scaladsl.server.Directives._
import cromiam.cromwell.CromwellClient
import cromwell.api.model._

trait WomtoolRouteSupport extends RequestSupport {
  val cromwellClient: CromwellClient

  val womtoolRoutes =
    path("womtool" / Segment / "describe") { _ =>
      post {
        extractStrictRequest { req =>
          complete {
            cromwellClient.forwardToCromwell(req).asHttpResponse
          }
        }
      }
    }

}
