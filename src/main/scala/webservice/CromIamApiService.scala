package webservice

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.server._
import akka.stream.ActorMaterializer
import spray.json._
import scala.concurrent.{ExecutionContextExecutor, Future}
import webservice.CromIamApiService._

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromiam"

  override def swaggerUiVersion = "2.1.1"
}

trait CromIamApiService extends Directives with SprayJsonSupport with DefaultJsonProtocol with RouteConcatenation {

  implicit val system: ActorSystem
  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  def cromwellInterface: String
  def cromwellPort: Int

  // Header that Akka HTTP adds to every request on receive.
  // We get an warning in logs if we don't strip it out before sending the request to cromwell
  // HTTP header ‘Timeout-Access: <function1>’ is not allowed in requests
  // See: https://github.com/akka/akka-http/issues/64
  val timeoutAccessHeader = "Timeout-Access"

  private[webservice] def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    val cromwellRequest = httpRequest
      .copy(uri = httpRequest.uri.withAuthority(cromwellInterface, cromwellPort))
      .withHeaders(httpRequest.headers.filterNot(header => header.name == timeoutAccessHeader))

    Http().singleRequest(cromwellRequest)
  }

  val workflowRoutes: Route = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute

  val engineRoutes: Route = statsRoute ~ versionRoute

  val allRoutes: Route = workflowRoutes ~ engineRoutes

  def statusRoute: Route =
    path("api" / "workflows" / Segment / Segment / "status") { (version, possibleWorkflowId) =>
      get {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def queryRoute: Route =
    path("api" / "workflows" / Segment / "query") { version =>
      parameterSeq { parameters =>
        get {
          extractRequest { req =>
            complete {
              forwardToCromwell(req)
            }
          }
        }
      }
    }

  def queryPostRoute: Route =
    path("api" / "workflows" / Segment / "query") { version =>
      (post & entity(as[Seq[Map[String, String]]])) { parameterMap =>
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def abortRoute: Route =
    path("api" / "workflows" / Segment / Segment / "abort") { (version, possibleWorkflowId) =>
      post {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def submitRoute: Route =
    path("api" / "workflows" / Segment) { version =>
      post {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def submitBatchRoute: Route =
    path("api" / "workflows" / Segment / "batch") { version =>
      post {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
        }
      }

  def workflowOutputsRoute: Route =
    path("api" / "workflows" / Segment / Segment / "outputs") { (version, possibleWorkflowId) =>
      get {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def workflowLogsRoute: Route =
    path("api" / "workflows" / Segment / Segment / "logs") { (version, possibleWorkflowId) =>
      get {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def metadataRoute: Route =
    path("api" / "workflows" / Segment / Segment / "metadata") { (version, possibleWorkflowId) =>
      get {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def timingRoute: Route =
    path("api" / "workflows" / Segment / Segment / "timing") { (version, possibleWorkflowId) =>
      get {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def statsRoute: Route =
    path("api" / "engine" / Segment / "stats") { version =>
      get {
        extractRequest { req =>
          complete {
            HttpResponse(status = Forbidden, entity = CromIamForbidden)
          }
        }
      }
    }

  def versionRoute: Route =
    path("api" / "engine" / Segment / "version") { version =>
      get {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }

  def backendRoute: Route =
    path("api" / "workflows" / Segment / "backends") { version =>
      get {
        extractRequest { req =>
          complete {
            forwardToCromwell(req)
          }
        }
      }
    }
}

object CromIamApiService {
  private[webservice] val CromIamForbidden = "CromIAM does not allow access to the /stats endpoint"
}
