package cromwell.webservice

import javax.ws.rs.Path

import akka.actor.{Actor, ActorRef, ActorRefFactory, Props}
import com.typesafe.config.Config
import com.wordnik.swagger.annotations._
import cromwell.engine.WorkflowId
import cromwell.engine.workflow.ValidateActor
import spray.http.StatusCodes
import spray.json._
import spray.routing.Directive.pimpApply
import spray.routing._

import scala.reflect.runtime.universe._
import scala.util.{Failure, Success, Try}

object SwaggerService {
  /*
    Because of the implicit arg requirement apply() doesn't work here, so falling back to the less
    idiomatic (but not unheard of) from().
   */
  def from(conf: Config)(implicit actorRefFactory: ActorRefFactory): SwaggerService = {
    new SwaggerService(conf.getConfig("swagger"))
  }
}

class SwaggerService(override val swaggerConfig: Config)
                    (implicit val actorRefFactory: ActorRefFactory)
  extends SwaggerConfigHttpService {
  override def apiTypes = Vector(typeOf[CromwellApiService])
}

object CromwellApiServiceActor {
  def props(workflowManagerActorRef: ActorRef, swaggerService: SwaggerService): Props = {
    Props(new CromwellApiServiceActor(workflowManagerActorRef, swaggerService))
  }
}

class CromwellApiServiceActor(val workflowManager: ActorRef, swaggerService: SwaggerService) extends Actor with CromwellApiService {
  implicit def executionContext = actorRefFactory.dispatcher
  def actorRefFactory = context

  def possibleRoutes = workflowRoutes ~ swaggerService.uiRoutes

  def receive = runRoute(possibleRoutes)
}

object CromwellApiService {
  /**
   * Used as swagger annotation constant below, this comma separated value lists the endpoint versions.
   * The last value is the latest version.
   */
  final val VersionAllowableValues = "v1"
}

trait CromwellApiService extends HttpService with PerRequestCreator {
  val workflowManager: ActorRef

  val workflowRoutes = queryRoute ~ workflowOutputsRoute ~ submitRoute ~ workflowStdoutStderrRoute ~ abortRoute ~
    callOutputsRoute ~ callStdoutStderrRoute ~ validateRoute

  def queryRoute =
    path("workflows" / Segment / Segment / "status") { (version, id) =>
      get {
        Try(WorkflowId.fromString(id)) match {
          case Success(workflowId) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowStatus(workflowId))
          case Failure(ex) =>
            complete(StatusCodes.BadRequest)
        }
      }
    }

  def abortRoute =
    path("workflows" / Segment / Segment / "abort") { (version, id) =>
      post {
        Try(WorkflowId.fromString(id)) match {
          case Success(workflowId) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowAbort(workflowId))
          case Failure(ex) =>
            complete(StatusCodes.BadRequest)
        }
      }
    }

  def submitRoute =
    path("workflows" / Segment) { version =>
      post {
        formFields("wdlSource", "workflowInputs") { (wdlSource, workflowInputs) =>

          val tryMap = Try(workflowInputs.parseJson)
          tryMap match {
            case Success(JsObject(json)) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowSubmit(wdlSource, workflowInputs, json))
            case Success(o) =>
              complete(StatusCodes.BadRequest, "Expecting JSON object as workflow inputs")
            case Failure(ex) =>
              complete(StatusCodes.BadRequest, "workflowInput JSON was malformed")
          }
        }
      }
    }

  def validateRoute =
    path("workflows" / Segment / "validate") { version =>
      post {
        formFields("wdlSource", "workflowInputs") { (wdlSource, workflowInputs) =>
          requestContext =>
            perRequest(
              requestContext,
              ValidateActor.props(wdlSource, workflowInputs),
              ValidateActor.ValidateWorkflow)
        }
      }
    }

  def workflowOutputsRoute =
    path("workflows" / Segment / Segment / "outputs") { (version, id) =>
      get {
        Try(WorkflowId.fromString(id)) match {
          case Success(workflowId) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowOutputs(workflowId))
          case Failure(ex) =>
            complete(StatusCodes.BadRequest)
        }
      }
    }

  def callOutputsRoute =
    path("workflows" / Segment / Segment / "outputs" / Segment) { (version, workflowId, callFqn) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(w) =>
          // This currently does not attempt to parse the call name for conformation to any pattern.
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.CallOutputs(w, callFqn))
        case Failure(_) =>
          complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
      }
    }

  def callStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs" / Segment) { (version, workflowId, callFqn) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(w) =>
          // This currently does not attempt to parse the call name for conformation to any pattern.
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.CallStdoutStderr(w, callFqn))
        case Failure(_) =>
          complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
      }
    }

  def workflowStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs") { (version, workflowId) =>
      Try(WorkflowId.fromString(workflowId)) match {
        case Success(w) =>
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowStdoutStderr(w))
        case Failure(_) =>
          complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
      }
    }
}
