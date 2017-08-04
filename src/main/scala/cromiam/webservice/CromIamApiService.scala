package cromiam.webservice

import akka.actor.ActorSystem
import akka.event.LoggingAdapter
import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.marshalling.ToResponseMarshallable
import akka.http.scaladsl.model.StatusCodes.{BadRequest, Forbidden, NotImplemented, Unauthorized, InternalServerError}
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.http.scaladsl.server._
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.ActorMaterializer
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cromiam.sam.SamClient
import cromiam.sam.SamClient.{SamDenialException, WorkflowAuthorizationRequest}
import cromiam.server.config.CromIamServerConfig
import cromiam.webservice.CromIamApiService._
import cromwell.api.model.CromwellStatus
import cromwell.api.model.CromwellStatusJsonSupport._
import spray.json._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContextExecutor, Future}

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromiam"
  override def swaggerUiVersion = "2.1.1"
}

trait CromIamApiService extends Directives with SprayJsonSupport with DefaultJsonProtocol with RouteConcatenation with SamClient {

  implicit val system: ActorSystem
  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  protected def configuration: CromIamServerConfig

  val logger: LoggingAdapter

  // Header that Akka HTTP adds to every request on receive.
  // We get an warning in logs if we don't strip it out before sending the request to cromwell
  // HTTP header ‘Timeout-Access: <function1>’ is not allowed in requests
  // See: https://github.com/akka/akka-http/issues/64
  val timeoutAccessHeader = "Timeout-Access"

  private[webservice] def forwardToCromwell(httpRequest: HttpRequest): Future[HttpResponse] = {
    val cromwellRequest = httpRequest
      .copy(uri = httpRequest.uri.withAuthority(configuration.cromwellConfig.interface, configuration.cromwellConfig.port).withScheme(configuration.cromwellConfig.scheme))
      .withHeaders(httpRequest.headers.filterNot(header => header.name == timeoutAccessHeader))
    Http().singleRequest(cromwellRequest)
  } recoverWith {
    case e => Future.failed(CromwellConnectionFailure(e))
  }

  private[webservice] def authorizeThenForwardToCromwell(authorization: Authorization, workflowIds: List[String], action: String, request: HttpRequest): Future[HttpResponse] = {
    def authForId(id: String): Future[Unit] = requestAuth(WorkflowAuthorizationRequest(authorization, id, action)) recoverWith {
      case SamDenialException => Future.failed(SamDenialException)
      case e => Future.failed(SamConnectionFailure("authorization", e))
    }

    (for {
      _ <- workflowIds traverse authForId
      resp <- forwardToCromwell(request)
    } yield resp) recover {
      case SamDenialException => SamDenialResponse
      case other => HttpResponse(status = InternalServerError, entity = s"CromIAM unexpected error: $other")
    }
  }

  private[webservice] def authorizeReadThenForwardToCromwell(authorization: Authorization, workflowId: List[String], request: HttpRequest): Future[HttpResponse] =
    authorizeThenForwardToCromwell(authorization, workflowId, "view", request)
  private[webservice] def authorizeAbortThenForwardToCromwell(authorization: Authorization, workflowId: String, request: HttpRequest): Future[HttpResponse] =
    authorizeThenForwardToCromwell(authorization, List(workflowId), "abort", request)

  private[webservice] def forwardSubmissionToCromwell(authorization: Authorization, request: HttpRequest, batch: Boolean): Future[HttpResponse] = {

    def extractWorkflowIds(response: HttpResponse): Future[List[String]] = {
      if (response.status.isSuccess) {
        if (batch) {
          Unmarshal(response).to[List[CromwellStatus]].map(_.map(_.id))
        } else {
          Unmarshal(response).to[CromwellStatus].map(s => List(s.id))
        }
      } else Future.successful(List.empty)
    }

    def registerWithSam(ids: List[String]): Future[Unit] = registerCreation(authorization, ids) recoverWith {
      case SamDenialException => Future.failed(SamDenialException)
      case e => Future.failed(SamConnectionFailure("new workflow registration", e))
    }

    (for {
      resp <- forwardToCromwell(request)
      workflowIds <- extractWorkflowIds(resp)
      _ <- registerWithSam(workflowIds)
    } yield resp) recover {
      case SamDenialException => SamDenialResponse
    }
  }

  val workflowRoutes: Route = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute ~ labelRoute ~ callCacheDiffRoute

  val engineRoutes: Route = statsRoute ~ versionRoute

  val allRoutes: Route = workflowRoutes ~ engineRoutes

  private def handleRequest(f: (Authorization, HttpRequest) => ToResponseMarshallable): Route = {
    headerValuePF { case a: Authorization => a } { authorization =>
      toStrictEntity(300.millis) {
        extractRequest { req =>
          complete {
            f.apply(authorization, req)
          }
        }
      }
    }
  }

  private def handleRequest(directive: Directive0)(f: (Authorization, HttpRequest) => ToResponseMarshallable): Route = {
    directive { handleRequest(f) }
  }

  def workflowGetRoute(urlSuffix: String): Route = workflowRoute(urlSuffix, get)

  def workflowRoute(urlSuffix: String, method: Directive0): Route = path("api" / "workflows" / Segment / Segment / urlSuffix) { (_, workflowId) =>
    handleRequest(method) { (userId, req) => authorizeReadThenForwardToCromwell(userId, List(workflowId), req) }
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
  def labelRoute: Route = workflowRoute("labels", patch)

  def backendRoute: Route = generalGetRoute("backends")

  def versionRoute: Route =  path("engine" / Segment / "version") { _ =>
    handleRequest(get) { (_, request) => forwardToCromwell(request) }
  }

  def statsRoute: Route = path("engine" / Segment / "stats") { _ =>
    handleRequest(get) { (_, _) => CromIamStatsForbidden }
  }

  def queryRoute: Route = path("api" / "workflows" / Segment / "query") { _ =>
    parameterSeq { _ =>
      handleRequest(get) { (_, _) => CromIamQueryNotImplemented }
    }
  }

  def queryPostRoute: Route = path("api" / "workflows" / Segment / "query") { _ =>
    (post & entity(as[Seq[Map[String, String]]])) { _ =>
      handleRequest { (_, _) => CromIamQueryNotImplemented }
    }
  }

  def callCacheDiffRoute: Route = path("api" / "workflows" / Segment / "callcaching" / "diff") { version =>
    parameterSeq { parameters =>
      handleRequest(get) { (user, req) =>
        val paramMap = parameters.toMap
        (paramMap.get("workflowA"), paramMap.get("workflowB")) match {
          case (Some(a), Some(b)) => authorizeReadThenForwardToCromwell(user, List(a, b), req)
          case _ => HttpResponse(status = BadRequest, entity = "Must supply both workflowA and workflowB to the /callcaching/diff endpoint")
        }
      }
    }
  }
}

object CromIamApiService {
  private[webservice] val SamDenialResponse = HttpResponse(status = Unauthorized, entity = "Access denied by Sam")
  private[webservice] val CromIamStatsForbidden = HttpResponse(status = Forbidden, entity = "CromIAM does not allow access to the /stats endpoint")
  private[webservice] val CromIamQueryNotImplemented = HttpResponse(status = NotImplemented, entity = "CromIAM does not currently support the /query endpoint")

  private[webservice] case class SamConnectionFailure(phase: String, f: Throwable) extends Exception(s"Unable to connect to Sam during $phase (${f.getMessage})", f)
  private[webservice] case class CromwellConnectionFailure(f: Throwable) extends Exception(s"Unable to connect to Cromwell (${f.getMessage})", f)
}
