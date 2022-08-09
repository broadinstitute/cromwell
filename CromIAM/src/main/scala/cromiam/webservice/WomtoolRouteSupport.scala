package cromiam.webservice

import akka.http.scaladsl.model.HttpResponse
import akka.http.scaladsl.server.Directives._
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient
import cromwell.api.model._

import scala.util.{Failure, Success}

//import scala.util.{Failure, Success}

trait WomtoolRouteSupport extends RequestSupport {
  // When this trait is mixed into `CromIamApiService` the value of `cromwellClient` is the reader (non-abort) address
  val cromwellClient: CromwellClient
  val samClient: SamClient

  val womtoolRoutes =
    path("api" / "womtool" / Segment / "describe") { _ =>
      post {
        extractUserAndRequest { (user, request) =>
          onComplete(samClient.isUserEnabledSam(user, request).value.unsafeToFuture()) {
            case Success(Left(httpResponse: HttpResponse)) => complete(httpResponse)
            case Success(Right(isEnabled: Boolean)) =>
              authorize(isEnabled) {
                complete {
                  cromwellClient.forwardToCromwell(request).asHttpResponse
                }
//                extractSubmission(user) { submission =>
//                  complete {
//                    forwardSubmissionToCromwell(
//                      user,
//                      submission.collection,
//                      request.withEntity(submission.entity)
//                    ).asHttpResponse
//                  }
//                }
              }
            case Failure(e) =>
              val message = s"Unable to look up submit whitelist for user ${user.userId}: ${e.getMessage}"
              throw new RuntimeException(message, e)

            // This endpoint requires authn which it gets for free from the proxy, does not care about authz
//            cromwellClient.forwardToCromwell(req).asHttpResponse
          }
        }
      }
    }

}
