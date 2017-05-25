package cromiam.webservice

import akka.actor.{ActorRef, ActorRefFactory, ActorSystem}
import akka.event.LoggingAdapter
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.marshalling.ToResponseMarshallable
import akka.http.scaladsl.model.StatusCodes.{Forbidden, NotImplemented, Unauthorized}
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.http.scaladsl.server._
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.ActorMaterializer
import cromiam.samiam.SamIamActor.{SamIamDenialException, WorkflowAuthorizationRequest}
import cromiam.samiam.{SamIamActor, SamIamClient}
import cromiam.server.config.CromIamServerConfig
import cromiam.webservice.CromIamApiService._
import cromwell.api.model.CromwellStatus
import cromwell.api.model.CromwellStatusJsonSupport._
import spray.json._

import scala.concurrent.{ExecutionContextExecutor, Future}
import scala.concurrent.duration._

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromiam"
  override def swaggerUiVersion = "2.1.1"
}

trait CromIamApiService extends Directives with SprayJsonSupport with DefaultJsonProtocol with RouteConcatenation with SamIamClient {

  implicit val system: ActorSystem
  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  protected def configuration: CromIamServerConfig
  protected lazy val UserIdHeader: String = configuration.cromIamConfig.userIdHeader

  override protected lazy val samIamActor: ActorRef = system.actorOf(SamIamActor.props(configuration.cromIamConfig.userIdHeader, configuration.cromIamConfig.allowedUsers), "SamIamActor")
  override protected lazy val actorRefFactory: ActorRefFactory = system

  val logger: LoggingAdapter

  // Header that Akka HTTP adds to every request on receive.
  // We get an warning in logs if we don't strip it out before sending the request to cromwell
  // HTTP header ‘Timeout-Access: <function1>’ is not allowed in requests
  // See: https://github.com/akka/akka-http/issues/64
  val timeoutAccessHeader = "Timeout-Access"

  private[webservice] def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    val cromwellRequest = httpRequest
      .copy(uri = httpRequest.uri.withAuthority(configuration.cromwellConfig.interface, configuration.cromwellConfig.port))
      .withHeaders(httpRequest.headers.filterNot(header => header.name == timeoutAccessHeader))

    Http().singleRequest(cromwellRequest)
  }

  private[webservice] def authorizeThenForwardToCromwell(user: Option[String], workflowId: String, action: String, request: HttpRequest): Future[HttpResponse] = {
    (for {
      _ <- requestAuth(WorkflowAuthorizationRequest(user.getOrElse("anon"), workflowId, action))
      resp <- forwardToCromwell(request)
    } yield resp) recover {
      case SamIamDenialException => HttpResponse(status = Unauthorized, entity = "")
    }
  }

  private[webservice] def authorizeReadThenForwardToCromwell(user: Option[String], workflowId: String, request: HttpRequest): Future[HttpResponse] =
    authorizeThenForwardToCromwell(user, workflowId, "read", request)
  private[webservice] def authorizeAbortThenForwardToCromwell(user: Option[String], workflowId: String, request: HttpRequest): Future[HttpResponse] =
    authorizeThenForwardToCromwell(user, workflowId, "abort", request)

  private[webservice] def forwardSubmissionToCromwell(user: Option[String], request: HttpRequest, batch: Boolean): Future[HttpResponse] = {

    def extractWorkflowIds(response: HttpResponse): Future[List[String]] = {
      if (response.status.isSuccess) {
        if (batch) {
          Unmarshal(response).to[List[CromwellStatus]].map(_.map(_.id))
        } else {
          Unmarshal(response).to[CromwellStatus].map(s => List(s.id))
        }
      } else Future.successful(List.empty)
    }

    val userNonOptional = user.getOrElse("anon")
    (for {
      resp <- forwardToCromwell(request)
      workflowIds <- extractWorkflowIds(resp)
      _ <- registerCreation(userNonOptional, workflowIds)
    } yield resp) recover {
      case SamIamDenialException => HttpResponse(status = Unauthorized, entity = "{}")
    }
  }

  val workflowRoutes: Route = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute

  val engineRoutes: Route = statsRoute ~ versionRoute

  val allRoutes: Route = workflowRoutes ~ engineRoutes

  private def handleRequest(f: (Option[String], HttpRequest) => ToResponseMarshallable): Route = {
    optionalHeaderValueByName(UserIdHeader) { userId =>
      toStrictEntity(300.millis) {
        extractRequest { req =>
          complete {
            f.apply(userId, req)
          }
        }
      }
    }
  }

  private def handleRequest(directive: Directive0)(f: (Option[String], HttpRequest) => ToResponseMarshallable): Route = {
    directive { handleRequest(f) }
  }

  def workflowGetRoute(urlSuffix: String): Route = path("api" / "workflows" / Segment / Segment / urlSuffix) { (_, workflowId) =>
    handleRequest(get) { (userId, req) => authorizeReadThenForwardToCromwell(userId, workflowId, req) }
  }
  def generalGetRoute(urlSuffix: String): Route = path("api" / "workflows" / Segment / urlSuffix) { _ =>
    handleRequest(get) { (_, req) => forwardToCromwell(req) }
  }

  def abortRoute: Route = path("api" / "workflows" / Segment / Segment / "abort") { (_, workflowId) =>
    handleRequest(post) { (userId, req) => authorizeAbortThenForwardToCromwell(userId, workflowId, req) }
  }

  def submitRoute: Route = path("api" / "workflows" / Segment) { _ =>
    handleRequest(post) { (userId, req) => forwardSubmissionToCromwell(userId, req, batch = false) }
  }
  def submitBatchRoute: Route = path("api" / "workflows" / Segment / "batch") { _ =>
    handleRequest(post) { (userId, req) => forwardSubmissionToCromwell(userId, req, batch = true) }
  }

  def workflowOutputsRoute: Route = workflowGetRoute("outputs")
  def workflowLogsRoute: Route = workflowGetRoute("logs")
  def metadataRoute: Route = workflowGetRoute("metadata")
  def timingRoute: Route = workflowGetRoute("metadata")
  def statusRoute: Route = workflowGetRoute("status")

  def versionRoute: Route = generalGetRoute("version")
  def backendRoute: Route = generalGetRoute("backends")

  def statsRoute: Route = path("api" / "engine" / Segment / "stats") { _ =>
    handleRequest(get) { (_, _) => CromIamStatsForbidden }
  }

  def queryRoute: Route = path("api" / "workflows" / Segment / "query") { version =>
    parameterSeq { parameters =>
      handleRequest(get) { (_, _) => CromIamQueryNotImplemented }
    }
  }

  def queryPostRoute: Route = path("api" / "workflows" / Segment / "query") { version =>
    (post & entity(as[Seq[Map[String, String]]])) { parameterMap =>
      handleRequest { (_, _) => CromIamQueryNotImplemented }
    }
  }
}

object CromIamApiService {
  private[webservice] val CromIamStatsForbidden = HttpResponse(status = Forbidden, entity = "CromIAM does not allow access to the /stats endpoint")
  private[webservice] val CromIamQueryNotImplemented = HttpResponse(status = NotImplemented, entity = "CromIAM does not currently support the /query endpoint")
}
