package cromiam.webservice

import akka.actor.ActorSystem
import akka.event.LoggingAdapter
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server._
import akka.http.scaladsl.server.Directives._
import akka.stream.ActorMaterializer
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient
import Collection.validateLabels
import cromiam.sam.SamClient
import cromiam.sam.SamClient.{CollectionAuthorizationRequest, SamConnectionFailure, SamDenialException, SamDenialResponse}
import cromiam.server.config.CromIamServerConfig
import cromiam.server.status.StatusService
import cromiam.webservice.CromIamApiService._

import scala.concurrent.{ExecutionContextExecutor, Future}

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromiam"
  override def swaggerUiVersion = "3.2.2" // TODO: Re-common-ize swagger out of cromwell's engine and reuse.
}

// NB: collection name *must* follow label value rules in cromwell. This needs to be documented somewhere. (although those restrictions are soon to die)
trait CromIamApiService extends RequestSupport
  with EngineRouteSupport
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
        log.error(e, "Request failed {}", e)
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


  val allRoutes: Route = handleExceptions(CromIamExceptionHandler) { workflowRoutes ~ engineRoutes }

  def abortRoute: Route = path("api" / "workflows" / Segment / Segment / Abort) { (_, workflowId) =>
    post {
      extractUserAndRequest { (user, req) =>
        logUserWorkflowAction(user, workflowId, Abort)
        complete {
          authorizeAbortThenForwardToCromwell(user, workflowId, req)
        }
      }
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
        extractUserAndRequest { (user, req) =>
          entity(as[String]) { labels =>
            logUserWorkflowAction(user, workflowId, Labels)
            validateLabels(Option(labels)) { _ => // Not using the labels, just using this to verify they didn't specify labels we don't want them to
              complete {
                authorizeReadThenForwardToCromwell(user, List(workflowId), req)
              }
            }
          }
        }
      }
    }
  }

  def backendRoute: Route = workflowGetRoute("backends")

  def callCacheDiffRoute: Route = path("api" / "workflows" / Segment / "callcaching" / "diff") { version =>
    get {
      extractUserAndRequest { (user, req) =>
        logUserAction(user, "call caching diff")
        parameterSeq { parameters =>
          val paramMap = parameters.toMap
          complete {
            (paramMap.get("workflowA"), paramMap.get("workflowB")) match {
              case (Some(a), Some(b)) => authorizeReadThenForwardToCromwell(user, List(a, b), req)
              case _ => HttpResponse(status = BadRequest, entity = "Must supply both workflowA and workflowB to the /callcaching/diff endpoint")
            }
          }
        }
      }
    }
  }

  /**
    * Base route for endpoints in the `workflows` space which do not take a workflow id as an argument
    */
  private def workflowGetRoute(urlSuffix: String): Route = path("api" / "workflows" / Segment / urlSuffix) { _ =>
    get {
      extractUserAndRequest { (user, req) =>
        logUserAction(user, urlSuffix)
        complete { cromwellClient.forwardToCromwell(req) }
      }
    }
  }

  /**
    * Base route for endpoints in the `workflows` space which take a workflow id as an argument
    */
  private def workflowGetRouteWithId(urlSuffix: String): Route = workflowRoute(urlSuffix, get)

  private def workflowRoute(urlSuffix: String, method: Directive0): Route = path("api" / "workflows" / Segment / Segment / urlSuffix) { (_, workflowId) =>
    method {
      extractUserAndRequest { (user, req) =>
        logUserWorkflowAction(user, workflowId, urlSuffix)
        complete { authorizeReadThenForwardToCromwell(user, List(workflowId), req) }
      }
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
          log.error(e, "Unable to connect to Sam {}", e)
          Future.failed(SamConnectionFailure("authorization", e))
      }
    }


    (for {
      rootWorkflowIds <- Future.sequence(workflowIds.map(id => cromwellClient.getRootWorkflow(id, user)))
      collections <- Future.sequence(rootWorkflowIds.map(id => cromwellClient.collectionForWorkflow(id, user))).map(_.distinct)
      _ <- collections traverse authForCollection
      resp <- cromwellClient.forwardToCromwell(request)
    } yield resp) recover {
      case SamDenialException => SamDenialResponse
      case other => HttpResponse(status = InternalServerError, entity = s"CromIAM unexpected error: $other")
    }
  }

  private def authorizeReadThenForwardToCromwell(user: User, workflowIds: List[String], request: HttpRequest): Future[HttpResponse] = {
    authorizeThenForwardToCromwell(user, workflowIds, "view", request)
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
