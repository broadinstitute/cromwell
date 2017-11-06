package cromiam.webservice

import akka.actor.ActorSystem
import akka.event.LoggingAdapter
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.marshalling.ToResponseMarshallable
import akka.http.scaladsl.model.StatusCodes.{BadRequest, Forbidden, InternalServerError, NotImplemented, Unauthorized}
import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.model.{HttpRequest, HttpResponse, StatusCodes}
import akka.http.scaladsl.server._
import akka.stream.ActorMaterializer
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient
import cromiam.sam.SamClient.{SamDenialException, WorkflowAuthorizationRequest}
import cromiam.server.config.CromIamServerConfig
import cromiam.server.status.StatusService
import cromiam.webservice.CromIamApiService._
import org.broadinstitute.dsde.workbench.util.health.StatusJsonSupport._
import spray.json._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContextExecutor, Future}

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromiam"
  override def swaggerUiVersion = "2.1.1"
}

trait CromIamApiService extends Directives with SprayJsonSupport with DefaultJsonProtocol with RouteConcatenation {

  implicit val system: ActorSystem
  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  protected def configuration: CromIamServerConfig

  lazy val cromwellClient = new CromwellClient(configuration.cromwellConfig.scheme, configuration.cromwellConfig.interface, configuration.cromwellConfig.port)
  lazy val samClient = new SamClient(configuration.samConfig.scheme, configuration.samConfig.interface, configuration.samConfig.port)

  val statusService: StatusService

  val logger: LoggingAdapter

  private[webservice] def authorizeThenForwardToCromwell(authorization: Authorization, workflowIds: List[String], action: String, request: HttpRequest): Future[HttpResponse] = {
    def authForId(id: String): Future[Unit] = samClient.requestAuth(WorkflowAuthorizationRequest(authorization, id, action)) recoverWith {
      case SamDenialException => Future.failed(SamDenialException)
      case e => Future.failed(SamConnectionFailure("authorization", e))
    }

    (for {
      _ <- workflowIds traverse authForId
      resp <- cromwellClient.forwardToCromwell(request)
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
    def registerWithSam(ids: List[String]): Future[Unit] = samClient.registerCreation(authorization, ids) recoverWith {
      case SamDenialException => Future.failed(SamDenialException)
      case e => Future.failed(SamConnectionFailure("new workflow registration", e))
    }

    (for {
      resp <- cromwellClient.forwardToCromwell(request)
      workflowIds <- cromwellClient.extractWorkflowIds(resp, batch)
      _ <- registerWithSam(workflowIds)
    } yield resp) recover {
      case SamDenialException => SamDenialResponse
    }
  }

  val workflowRoutes: Route = queryRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~ submitBatchRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute ~ labelRoute ~ callCacheDiffRoute

  val engineRoutes: Route = statsRoute ~ versionRoute ~ engineStatusRoute

  val allRoutes: Route = workflowRoutes ~ engineRoutes

  private def handleRequestWithAuthn(directive: Directive0)(f: (Authorization, HttpRequest) => ToResponseMarshallable): Route = {
    directive {
      headerValuePF { case a: Authorization => a } { authorization =>
        toStrictEntity(Timeout) {
          extractRequest { req =>
            complete {
              f.apply(authorization, req)
            }
          }
        }
      }
    }
  }

  private def handlePublicRequest(directive: Directive0)(f: (HttpRequest) => ToResponseMarshallable): Route = {
    directive {
      toStrictEntity(Timeout) {
        extractRequest { req =>
          complete {
            f.apply(req)
          }
        }
      }
    }
  }

  /**
    * Base route for endpoints in the `workflows` space which do not take a workflow id as an argument
    */
  def workflowGetRoute(urlSuffix: String): Route = path("api" / "workflows" / Segment / urlSuffix) { _ =>
    handleRequestWithAuthn(get) { (_, req) => cromwellClient.forwardToCromwell(req) }
  }

  /**
    * Base route for endpoints in the `workflows` space which take a workflow id as an argument
    */
  def workflowGetRouteWithId(urlSuffix: String): Route = workflowRoute(urlSuffix, get)

  def workflowRoute(urlSuffix: String, method: Directive0): Route = path("api" / "workflows" / Segment / Segment / urlSuffix) { (_, workflowId) =>
    handleRequestWithAuthn(method) { (userId, req) => authorizeReadThenForwardToCromwell(userId, List(workflowId), req) }
  }

  def abortRoute: Route = path("api" / "workflows" / Segment / Segment / "abort") { (_, workflowId) =>
    handleRequestWithAuthn(post) { (userId, req) => authorizeAbortThenForwardToCromwell(userId, workflowId, req) }
  }

  def submitRoute: Route = path("api" / "workflows" / Segment) { _ =>
    handleRequestWithAuthn(post) { (userId, req) => forwardSubmissionToCromwell(userId, req, batch = false) }
  }
  def submitBatchRoute: Route = path("api" / "workflows" / Segment / "batch") { _ =>
    handleRequestWithAuthn(post) { (userId, req) => forwardSubmissionToCromwell(userId, req, batch = true) }
  }

  def workflowOutputsRoute: Route = workflowGetRouteWithId("outputs")
  def workflowLogsRoute: Route = workflowGetRouteWithId("logs")
  def metadataRoute: Route = workflowGetRouteWithId("metadata")
  def timingRoute: Route = workflowGetRouteWithId("metadata")
  def statusRoute: Route = workflowGetRouteWithId("status")
  def labelRoute: Route = workflowRoute("labels", patch)

  def backendRoute: Route = workflowGetRoute("backends")

  def versionRoute: Route = path("engine" / Segment / "version") { _ => handlePublicRequest(get) { req => cromwellClient.forwardToCromwell(req) } }

  def engineStatusRoute: Route = path("engine" / Segment / "status") { _ =>
    get {
      complete(statusService.status().map { statusResponse =>
        val httpStatus = if (statusResponse.ok) StatusCodes.OK else StatusCodes.InternalServerError
        (httpStatus, statusResponse)
      })
    }
  }

  def statsRoute: Route = path("engine" / Segment / "stats") { _ => handlePublicRequest(get) { _ => CromIamStatsForbidden } }

  def queryRoute: Route = path("api" / "workflows" / Segment / "query") { _ =>
    parameterSeq { _ =>  handleRequestWithAuthn(get) { (_, _) =>  CromIamQueryNotImplemented } }
  }

  def queryPostRoute: Route = path("api" / "workflows" / Segment / "query") { _ =>
    (post & entity(as[Seq[Map[String, String]]])) { _ =>
      handleRequestWithAuthn(post) { (_, _) => CromIamQueryNotImplemented }
    }
  }

  def callCacheDiffRoute: Route = path("api" / "workflows" / Segment / "callcaching" / "diff") { version =>
    parameterSeq { parameters =>
      handleRequestWithAuthn(get) { (user, req) =>
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

  def Timeout = 300.millis
}
