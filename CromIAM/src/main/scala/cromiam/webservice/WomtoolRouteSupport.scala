package cromiam.webservice

import akka.http.scaladsl.server.Directives._
import cromiam.cromwell.CromwellClient
import cromwell.api.model._

trait WomtoolRouteSupport extends RequestSupport {
  // When this trait is mixed into `CromIamApiService` the value of `cromwellClient` is the reader (non-abort) address
  val cromwellClient: CromwellClient

  val womtoolRoutes =
    path("api" / "womtool" / Segment / "describe") { _ =>
      post {
        extractStrictRequest { req =>
          complete {
            // This endpoint requires authn which it gets for free from the proxy, does not care about authz
            cromwellClient.forwardToCromwell(req).asHttpResponse
          }
        }
      }
    }

}
