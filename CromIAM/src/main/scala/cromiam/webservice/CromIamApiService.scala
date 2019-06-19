package cromiam.webservice

import akka.actor.{ActorRef, ActorSystem}
import akka.event.LoggingAdapter
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import akka.stream.ActorMaterializer
import cats.effect.IO
import cats.instances.list._
import cats.syntax.traverse._
import com.typesafe.config.Config
import cromiam.auth.Collection.validateLabels
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient
import cromiam.instrumentation.CromIamInstrumentation
import cromiam.sam.SamClient
import cromiam.sam.SamClient.{CollectionAuthorizationRequest, SamConnectionFailure, SamDenialException, SamDenialResponse}
import cromiam.server.config.CromIamServerConfig
import cromiam.server.status.StatusService
import cromiam.webservice.CromIamApiService._
import cromwell.api.model._
import cromwell.services.ServiceRegistryActor

import scala.concurrent.ExecutionContextExecutor

trait SwaggerService extends SwaggerUiResourceHttpService {
  override def swaggerServiceName = "cromiam"
  override def swaggerUiVersion = "3.2.2" // TODO: Re-common-ize swagger out of cromwell's engine and reuse.
}

// NB: collection name *must* follow label value rules in cromwell. This needs to be documented somewhere. (although those restrictions are soon to die)
trait CromIamApiService extends RequestSupport
  with EngineRouteSupport
  with SubmissionSupport
  with QuerySupport
  with WomtoolRouteSupport
  with CromIamInstrumentation {

  implicit val system: ActorSystem
  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  protected def rootConfig: Config
  protected def configuration: CromIamServerConfig

  override lazy val serviceRegistryActor: ActorRef = system.actorOf(ServiceRegistryActor.props(rootConfig), "ServiceRegistryActor")

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
    log,
    serviceRegistryActor)

  lazy val cromwellAbortClient = new CromwellClient(configuration.cromwellAbortConfig.scheme,
    configuration.cromwellAbortConfig.interface,
    configuration.cromwellAbortConfig.port,
    log,
    serviceRegistryActor)

  lazy val samClient = new SamClient(
    configuration.samConfig.http.scheme,
    configuration.samConfig.http.interface,
    configuration.samConfig.http.port,
    configuration.samConfig.checkSubmitWhitelist,
    log,
    serviceRegistryActor)

  val statusService: StatusService

  val workflowRoutes: Route = queryGetRoute ~ queryPostRoute ~ workflowOutputsRoute ~ submitRoute ~
    workflowLogsRoute ~ abortRoute ~ metadataRoute ~ timingRoute ~ statusRoute ~ backendRoute ~ labelPatchRoute ~
    callCacheDiffRoute ~ labelGetRoute ~ releaseHoldRoute


  val allRoutes: Route = handleExceptions(CromIamExceptionHandler) { workflowRoutes ~ engineRoutes ~ womtoolRoutes }

  def abortRoute: Route = path("api" / "workflows" / Segment / Segment / Abort) { (_, workflowId) =>
    post {
      extractUserAndRequest { (user, req) =>
        logUserWorkflowAction(user, workflowId, Abort)
        complete {
          authorizeAbortThenForwardToCromwell(user, workflowId, req).asHttpResponse
        }
      }
    }
  }

  //noinspection MutatorLikeMethodIsParameterless
  def releaseHoldRoute: Route =  path("api" / "workflows" / Segment / Segment / ReleaseHold) { (_, workflowId) =>
    post {
      extractUserAndRequest { (user, req) =>
        logUserWorkflowAction(user, workflowId, ReleaseHold)
        complete {
          authorizeUpdateThenForwardToCromwell(user, workflowId, req).asHttpResponse
        }
      }
    }
  }

  def workflowOutputsRoute: Route = workflowGetRouteWithId("outputs")
  def workflowLogsRoute: Route = workflowGetRouteWithId("logs")
  def metadataRoute: Route = workflowGetRouteWithId("metadata")
  def timingRoute: Route = workflowGetRouteWithId("timing")
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
                authorizeUpdateThenForwardToCromwell(user, workflowId, req).asHttpResponse
              }
            }
          }
        }
      }
    }
  }

  def backendRoute: Route = workflowGetRoute("backends")

  def callCacheDiffRoute: Route = path("api" / "workflows" / Segment / "callcaching" / "diff") { _ =>
    get {
      extractUserAndRequest { (user, req) =>
        logUserAction(user, "call caching diff")
        parameterSeq { parameters =>
          val paramMap = parameters.toMap
          complete {
            (paramMap.get("workflowA"), paramMap.get("workflowB")) match {
              case (Some(a), Some(b)) => authorizeReadThenForwardToCromwell(user, List(a, b), req).asHttpResponse
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
        complete {
          cromwellClient.forwardToCromwell(req).asHttpResponse
        }
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
        complete {
          authorizeReadThenForwardToCromwell(user, List(workflowId), req).asHttpResponse
        }
      }
    }
  }

  /**
    * Authorization-related Cromwell interactions run through the `cromwellAuthClient`, while the actual Cromwell
    * request runs through `cromwellRequestClient`. These can be the same Cromwell instance but they can also be
    * different if the request needs to go to a specific Cromwell instance such as an abort server.
    */
  private def authorizeThenForwardToCromwell(user: User,
                                             workflowIds: List[String],
                                             action: String,
                                             request: HttpRequest,
                                             cromwellAuthClient: CromwellClient,
                                             cromwellRequestClient: CromwellClient):
  FailureResponseOrT[HttpResponse] = {
    def authForCollection(collection: Collection): FailureResponseOrT[Unit] = {
      samClient.requestAuth(CollectionAuthorizationRequest(user, collection, action), request) mapErrorWith {
        case e: SamDenialException => IO.raiseError(e)
        case e =>
          log.error(e, "Unable to connect to Sam {}", e)
          IO.raiseError(SamConnectionFailure("authorization", e))
      }
    }

    val cromwellResponseT = for {
      rootWorkflowIds <- workflowIds.traverse(cromwellAuthClient.getRootWorkflow(_, user, request))
      collections <- rootWorkflowIds
        .traverse(cromwellAuthClient.collectionForWorkflow(_, user, request))
        .map(_.distinct)
      _ <- collections traverse authForCollection
      resp <- cromwellRequestClient.forwardToCromwell(request)
    } yield resp

    FailureResponseOrT(
      cromwellResponseT.value handleErrorWith {
        case _: SamDenialException => IO.pure(Left(SamDenialResponse))
        case other =>
          IO.pure(Left(HttpResponse(status = InternalServerError, entity = s"CromIAM unexpected error: $other")))
      }
    )
  }

  private def authorizeReadThenForwardToCromwell(user: User,
                                                 workflowIds: List[String],
                                                 request: HttpRequest
                                                ): FailureResponseOrT[HttpResponse] = {
    authorizeThenForwardToCromwell(
      user = user,
      workflowIds = workflowIds,
      action = "view",
      request = request,
      cromwellAuthClient = cromwellClient,
      cromwellRequestClient = cromwellClient)
  }

  private def authorizeUpdateThenForwardToCromwell(user: User,
                                                   workflowId: String,
                                                   request: HttpRequest
                                                  ): FailureResponseOrT[HttpResponse] = {
    authorizeThenForwardToCromwell(
      user = user,
      workflowIds = List(workflowId),
      action = "update",
      request = request,
      cromwellAuthClient = cromwellClient,
      cromwellRequestClient = cromwellClient)
  }

  private def authorizeAbortThenForwardToCromwell(user: User,
                                                  workflowId: String,
                                                  request: HttpRequest
                                                 ): FailureResponseOrT[HttpResponse] = {
    // Do all the authing for the abort with "this" cromwell instance (cromwellClient), but the actual abort command
    // must go to the dedicated abort server (cromwellAbortClient).
    authorizeThenForwardToCromwell(
      user = user,
      workflowIds = List(workflowId),
      action = "abort",
      request = request,
      cromwellAuthClient = cromwellClient,
      cromwellRequestClient = cromwellAbortClient
    )
  }

  private def logUserAction(user: User, action: String) = log.info("User " + user.userId + " requesting " + action)

  private def logUserWorkflowAction(user: User, wfId: String, action: String) = {
    log.info("User " + user.userId + " requesting " + action + " with " + wfId)
  }
}

object CromIamApiService {
  val Abort = "abort"
  val Labels = "labels"
  val ReleaseHold = "releaseHold"
}
