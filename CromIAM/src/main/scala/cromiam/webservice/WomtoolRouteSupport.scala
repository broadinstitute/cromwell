package cromiam.webservice

import akka.http.scaladsl.server.Directives._
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient

trait WomtoolRouteSupport extends RequestSupport {
  // When this trait is mixed into `CromIamApiService` the value of `cromwellClient` is the reader (non-abort) address
  val cromwellClient: CromwellClient
  val samClient: SamClient

  val womtoolRoutes =
    path("api" / "womtool" / Segment / "describe") { _ =>
      post {
        extractUserAndStrictRequest { (user, req) =>
          forwardIfUserEnabled(user, req, cromwellClient, samClient)
        }
      }
    }

}
