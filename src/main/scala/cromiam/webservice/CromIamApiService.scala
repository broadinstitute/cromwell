package cromiam.webservice

import akka.actor.ActorSystem
import akka.event.LoggingAdapter
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import akka.stream.ActorMaterializer
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient
import Collection.validateLabels
import cromiam.sam.SamClient
import cromiam.sam.SamClient.{CollectionAuthorizationRequest, SamConnectionFailure, SamDenialException, SamDenialResponse}
import cromiam.webservice.RequestDirectives._
import cromiam.server.config.CromIamServerConfig
import cromiam.server.status.StatusService
import cromiam.webservice.CromIamApiService._
import org.broadinstitute.dsde.workbench.util.health.StatusJsonSupport._
import spray.json._

import scala.concurrent.{ExecutionContextExecutor, Future}

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromiam"
  override def swaggerUiVersion = "2.1.1"
}

// NB: collection name *must* follow label value rules in cromwell. This needs to be documented somewhere. (although those restrictions are soon to die)
trait CromIamApiService extends Directives
  with SprayJsonSupport
  with DefaultJsonProtocol
  with RouteConcatenation
  with SubmissionSupport
  with QuerySupport {

  implicit val system: ActorSystem
  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  protected def configuration: CromIamServerConfig

  val log: LoggingAdapter

  val CromIamExceptionHandler: ExceptionHandler = {
  ExceptionHandler {
      case e: Exception =>
        log.error("Request failed", e)
        complete(HttpResponse(InternalServerError, entity = e.getMessage)) // FIXME: use workbench-model ErrorReport
    }
  }

  lazy val cromwellClient = new CromwellClient(configuration.cromwellConfig.scheme,
    configuration.cromwellConfig.interface,
    configuration.cromwellConfig.port,
    log)

  lazy val samClient = new SamClient(configuration.samConfig.scheme, configuration.samConfig.interface, configuration.samConfig.port, log)

  val statusService: StatusService

  val workflowRoutes: Route = queryGetRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute ~ labelPatchRoute ~
    callCacheDiffRoute ~ labelGetRoute

  val engineRoutes: Route = statsRoute ~ versionRoute ~ engineStatusRoute

  val allRoutes: Route = handleExceptions(CromIamExceptionHandler) { workflowRoutes ~ engineRoutes }

  def abortRoute: Route = path("api" / "workflows" / Segment / Segment / Abort) { (_, workflowId) =>
    handleRequestWithAuthn(post) { (user, req) =>
      logUserWorkflowAction(user, workflowId, Abort)
      authorizeAbortThenForwardToCromwell(user, workflowId, req)
    }
  }

  def workflowOutputsRoute: Route = workflowGetRouteWithId("outputs")
  def workflowLogsRoute: Route = workflowGetRouteWithId("logs")
  def metadataRoute: Route = workflowGetRouteWithId("metadata")
  def timingRoute: Route = workflowGetRouteWithId("metadata")
  def statusRoute: Route = workflowGetRouteWithId("status")
  def labelGetRoute: Route = workflowGetRouteWithId(Labels)

  def labelPatchRoute: Route = {
    path("api" / "workflows" / Segment / Segment / Labels) { (_, workflowId) =>
      patch {
        entity(as[String]) { labels =>
          validateLabels(Option(labels)) { _ => // Not using the labels, just using this to verify they didn't specify labels we don't want them to
            handleRequestWithAuthn(patch) { (user, req) =>
              logUserWorkflowAction(user, workflowId, Labels)
              authorizeReadThenForwardToCromwell(user, List(workflowId), req)
            }
          }
        }
      }
    }
  }

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

  def callCacheDiffRoute: Route = path("api" / "workflows" / Segment / "callcaching" / "diff") { version =>
    parameterSeq { parameters =>
      handleRequestWithAuthn(get) { (user, req) =>
        logUserAction(user, "call caching diff")
        val paramMap = parameters.toMap
        (paramMap.get("workflowA"), paramMap.get("workflowB")) match {
          case (Some(a), Some(b)) => authorizeReadThenForwardToCromwell(user, List(a, b), req)
          case _ => HttpResponse(status = BadRequest, entity = "Must supply both workflowA and workflowB to the /callcaching/diff endpoint")
        }
      }
    }
  }

  /**
    * Base route for endpoints in the `workflows` space which do not take a workflow id as an argument
    */
  private def workflowGetRoute(urlSuffix: String): Route = path("api" / "workflows" / Segment / urlSuffix) { _ =>
    handleRequestWithAuthn(get) { (user, req) =>
      logUserAction(user, urlSuffix)
      cromwellClient.forwardToCromwell(req)
    }
  }

  /**
    * Base route for endpoints in the `workflows` space which take a workflow id as an argument
    */
  private def workflowGetRouteWithId(urlSuffix: String): Route = workflowRoute(urlSuffix, get)

  private def workflowRoute(urlSuffix: String, method: Directive0): Route = path("api" / "workflows" / Segment / Segment / urlSuffix) { (_, workflowId) =>
    handleRequestWithAuthn(method) { (user, req) =>
      logUserWorkflowAction(user, workflowId, urlSuffix)
      authorizeReadThenForwardToCromwell(user, List(workflowId), req)
    }
  }

  private def authorizeThenForwardToCromwell(user: User,
                                             workflowIds: List[String],
                                             action: String,
                                             request: HttpRequest): Future[HttpResponse] = {
    def authForCollection(collection: Collection): Future[Unit] = {
      samClient.requestAuth(CollectionAuthorizationRequest(user, collection, action)) recoverWith {
        case SamDenialException => Future.failed(SamDenialException)
        case e =>
          log.error("Unable to connect to Sam", e)
          Future.failed(SamConnectionFailure("authorization", e))
      }
    }

    (for {
      collections <- Future.sequence(workflowIds.map(id => cromwellClient.collectionForWorkflow(id, user))).map(_.distinct)
      _ <- collections traverse authForCollection
      resp <- cromwellClient.forwardToCromwell(request)
    } yield resp) recover {
      case SamDenialException => SamDenialResponse
      case other => HttpResponse(status = InternalServerError, entity = s"CromIAM unexpected error: $other")
    }
  }

  // FIXME: Refactor target - why the heck does one take a List[String] and one take a String and List-ify it? Leaving as-is for now
  private def authorizeReadThenForwardToCromwell(user: User, workflowId: List[String], request: HttpRequest): Future[HttpResponse] = {
    authorizeThenForwardToCromwell(user, workflowId, "view", request)
  }

  private def authorizeAbortThenForwardToCromwell(user: User, workflowId: String, request: HttpRequest): Future[HttpResponse] = {
    authorizeThenForwardToCromwell(user, List(workflowId), "abort", request)
  }

  private def logUserAction(user: User, action: String) = log.info("User " + user.userId + " requesting " + action)

  private def logUserWorkflowAction(user: User, wfId: String, action: String) = {
    log.info("User " + user.userId + " requesting " + action + " with " + wfId)
  }
}

object CromIamApiService {
  private[webservice] val CromIamStatsForbidden = HttpResponse(status = Forbidden, entity = "CromIAM does not allow access to the /stats endpoint")

  val Abort = "abort"
  val Labels = "labels"
}
