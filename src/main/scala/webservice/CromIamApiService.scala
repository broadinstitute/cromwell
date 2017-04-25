package webservice

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.model.HttpResponse
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.server._
import spray.json._

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromiam"

  override def swaggerUiVersion = "2.1.1"
}

trait CromIamApiService extends Directives with SprayJsonSupport with DefaultJsonProtocol with RouteConcatenation {

  def returnInternalServerError(msg: String)  = HttpResponse(InternalServerError, entity = msg)

  val workflowRoutes = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute

  val engineRoutes = statsRoute ~ versionRoute

  val allRoutes = workflowRoutes ~ engineRoutes

  def statusRoute =
    path("api" / "workflows" / Segment / Segment / "status") { (version, possibleWorkflowId) =>
      get {
        complete {
          returnInternalServerError("workflow status")
        }
      }
    }

  def queryRoute =
    path("api" / "workflows" / Segment / "query") { version =>
      parameterSeq { parameters =>
        get {
          complete {
            returnInternalServerError("query get")
          }
        }
      }
    }

  def queryPostRoute =
    path("api" / "workflows" / Segment / "query") { version =>
      (post & entity(as[Seq[Map[String, String]]])) { parameterMap =>
          complete {
            returnInternalServerError("workflow query post")
          }
      }
    }

  def abortRoute =
    path("api" / "workflows" / Segment / Segment / "abort") { (version, possibleWorkflowId) =>
      post {
        complete {
          returnInternalServerError("workflow abort")
        }
      }
    }

  def submitRoute =
    path("api" / "workflows" / Segment) { version =>
      post {
        complete{
          returnInternalServerError("submit workflow")
        }
      }
    }

  def submitBatchRoute =
    path("api" / "workflows" / Segment / "batch") { version =>
      post {
        complete {
          returnInternalServerError("batch submit workflow")
          }
        }
      }

  def workflowOutputsRoute =
    path("api" / "workflows" / Segment / Segment / "outputs") { (version, possibleWorkflowId) =>
      get {
        complete {
          returnInternalServerError("workflow outputs")
        }
      }
    }

  def workflowLogsRoute =
    path("api" / "workflows" / Segment / Segment / "logs") { (version, possibleWorkflowId) =>
      get {
        complete {
          returnInternalServerError("workflow logs")
        }
      }
    }

  def metadataRoute =
    path("api" / "workflows" / Segment / Segment / "metadata") { (version, possibleWorkflowId) =>
      get {
        complete {
          returnInternalServerError("workflow metdata")
        }
      }
    }

  def timingRoute =
    path("api" / "workflows" / Segment / Segment / "timing") { (version, possibleWorkflowId) =>
      get {
        complete {
          returnInternalServerError("workflow timing")
        }
      }
    }

  def statsRoute =
    path("api" / "engine" / Segment / "stats") { version =>
      get {
        complete {
          returnInternalServerError("engine stats")
        }
      }
    }

  def versionRoute =
    path("api" / "engine" / Segment / "version") { version =>
      get {
        complete {
          returnInternalServerError("engine version")
        }
      }
    }

  def backendRoute =
    path("api" / "workflows" / Segment / "backends") { version =>
      get {
        complete {
          returnInternalServerError("engine backends")
        }
      }
    }

}
