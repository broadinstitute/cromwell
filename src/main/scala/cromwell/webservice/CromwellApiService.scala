package cromwell.webservice

import java.util.UUID
import javax.ws.rs.Path

import akka.actor.{Actor, ActorRef, ActorRefFactory, Props}
import com.typesafe.config.Config
import com.wordnik.swagger.annotations._
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

  def possibleRoutes = workflowRoutes ~ swaggerService.routes ~ swaggerService.uiRoutes

  def receive = runRoute(possibleRoutes)
}

object CromwellApiService {
  /**
   * Used as swagger annotation constant below, this comma separated value lists the endpoint versions.
   * The last value is the latest version.
   */
  final val VersionAllowableValues = "v1"
}

/*
 NOTE: Please run `sbt 'run server'` and regression test by hand the http://localhost:8000/swagger UI.

 "The" swagger implementation consists of a number of components: the spec, the core, the ui, and the spray wrappers.
 spray-swagger uses swagger-core to generate swagger-spec compliant JSON, that swagger-ui then renders. All components
 currently support the baseline of spec 1.2, but there are a number of bugs/issues related to how they interact with
 each other.

 Examples:
  Setting a `response` with class type of a scala `Map` seemed to cause the entire UI not to render.
  The spec says a `dataType` of `"File"` must be uppercase for form uploads.
  etc.

 There is work in progress to upgrade the stack to swagger-spec 2.0, but it appears to be tied to akka-http not spray.
 */

@Api(value = "/workflows", description = "Workflow Services", produces = "application/json", position = 1)
trait CromwellApiService extends HttpService with PerRequestCreator {
  val workflowManager: ActorRef

  val workflowRoutes = queryRoute ~ workflowOutputsRoute ~ submitRoute ~ callOutputsRoute ~ callStdoutStderrRoute ~ abortRoute

  @Path("/{version}/{id}/status")
  @ApiOperation(
    value = "Query for workflow status based on workflow id",
    nickname = "status",
    httpMethod = "GET",
    produces = "application/json"
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "version", required = true, dataType = "string", paramType = "path", value = "API Version", allowableValues = CromwellApiService.VersionAllowableValues),
    new ApiImplicitParam(name = "id", required = true, dataType = "string", paramType = "path", value = "Workflow Identifier")
  ))
  @ApiResponses(Array(
    new ApiResponse(code = 200, message = "Successful Request", response = classOf[WorkflowStatusResponse]),
    new ApiResponse(code = 404, message = "Workflow ID Not Found"),
    new ApiResponse(code = 400, message = "Malformed Workflow ID"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def queryRoute =
    path("workflows" / Segment / Segment / "status") { (version, id) =>
      get {
        Try(UUID.fromString(id)) match {
          case Success(workflowId) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowStatus(workflowId))
          case Failure(ex) =>
            complete(StatusCodes.BadRequest)
        }
      }
    }

  @Path("/{version}/{id}/abort")
  @ApiOperation(
    value = "Abort a workflow based on workflow id.",
    nickname = "abort",
    httpMethod = "POST",
    produces = "application/json"
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "version", required = true, dataType = "string", paramType = "path", value = "API Version", allowableValues = CromwellApiService.VersionAllowableValues),
    new ApiImplicitParam(name = "id", required = true, dataType = "string", paramType = "path", value = "workflow identifier")
  ))
  @ApiResponses(Array(
    new ApiResponse(code = 200, message = "Successful Request", response = classOf[WorkflowAbortResponse]),
    new ApiResponse(code = 404, message = "Workflow ID Not Found"),
    new ApiResponse(code = 403, message = "Workflow in terminal status"),
    new ApiResponse(code = 400, message = "Malformed Workflow ID"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def abortRoute =
    path("workflows" / Segment / Segment / "abort") { (version, id) =>
      post {
        Try(UUID.fromString(id)) match {
          case Success(workflowId) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowAbort(workflowId))
          case Failure(ex) =>
            complete(StatusCodes.BadRequest)
        }
      }
    }

  @ApiOperation(
    value = "Submit a new workflow for execution",
    nickname = "submit",
    httpMethod = "POST",
    consumes = "multipart/form-data",
    produces = "application/json"
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "version", required = true, dataType = "string", paramType = "path", value = "API Version", allowableValues = CromwellApiService.VersionAllowableValues),
    new ApiImplicitParam(name = "wdlSource", required = true, dataType = "File", paramType = "form", value = "WDL Source"),
    new ApiImplicitParam(name = "workflowInputs", required = true, dataType = "File", paramType = "form", value = "WDL JSON")
  ))
  @ApiResponses(Array(
    new ApiResponse(code = 201, message = "Successful Request", response = classOf[WorkflowStatusResponse]),
    new ApiResponse(code = 400, message = "Malformed Input"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
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

  @Path("/{version}/{id}/outputs")
  @ApiOperation(
    value = "Query for workflow outputs based on workflow id.",
    nickname = "workflow-outputs",
    httpMethod = "GET",
    produces = "application/json"
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "version", required = true, dataType = "string", paramType = "path", value = "API Version", allowableValues = CromwellApiService.VersionAllowableValues),
    new ApiImplicitParam(name = "id", required = true, dataType = "string", paramType = "path", value = "Workflow Identifier")
  ))
  @ApiResponses(Array(
    // NOTE: Uncommenting the `response = classOf[WorkflowOutputs]` below with the current swagger stack causes the
    // entire "Workflow" section to disappear from the swagger ui.
    // Non-basic types, including Maps, seem to be very broken throughout the swagger stack:
    // ex: https://github.com/swagger-api/swagger-core/issues/702
    new ApiResponse(code = 200, message = "Successful Request"),
    // breaks - new ApiResponse(code = 200, message = "Successful Request", response = classOf[WorkflowOutputs]),
    // ugly   - new ApiResponse(code = 200, message = "Successful Request", response = classOf[java.util.Map[String, WdlValue]]),
    new ApiResponse(code = 404, message = "Workflow ID Not Found"),
    new ApiResponse(code = 400, message = "Malformed Workflow ID"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def workflowOutputsRoute =
    path("workflows" / Segment / Segment / "outputs") { (version, id) =>
      get {
        Try(UUID.fromString(id)) match {
          case Success(workflowId) =>
            requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowOutputs(workflowId))
          case Failure(ex) =>
            complete(StatusCodes.BadRequest)
        }
      }
    }

  @Path("/{version}/{workflowId}/outputs/{callFqn}")
  @ApiOperation(
    value = "Query for call outputs based on workflow id and call fully qualified name (e.g. my_workflow.my_call).",
    nickname = "call-outputs",
    httpMethod = "GET",
    produces = "application/json"
    //,response = classOf[CallOutputs]
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "version", required = true, dataType = "string", paramType = "path", value = "API Version", allowableValues = CromwellApiService.VersionAllowableValues),
    new ApiImplicitParam(name = "workflowId", required = true, dataType = "string", paramType = "path", value = "Workflow ID"),
    new ApiImplicitParam(name = "callFqn", required = true, dataType = "string", paramType = "path", value = "Call fully qualified name")
  ))
  @ApiResponses(Array(
    // See the note above regarding Swagger brokenness for why there is no response classOf here.
    new ApiResponse(code = 200, message = "Successful Request"),
    new ApiResponse(code = 404, message = "Workflow ID Not Found"),
    new ApiResponse(code = 404, message = "Call Fully Qualified Name Not Found"),
    new ApiResponse(code = 400, message = "Malformed Workflow ID"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def callOutputsRoute =
    path("workflows" / Segment / Segment / "outputs" / Segment) { (version, workflowId, callFqn) =>
      Try(UUID.fromString(workflowId)) match {
        case Success(w) =>
          // This currently does not attempt to parse the call name for conformation to any pattern.
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.CallOutputs(w, callFqn))
        case Failure(ex) =>
          complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
      }
    }

  @Path("/{version}/{workflowId}/stdout/{callFqn}")
  @ApiOperation(
    value = "Query for the standard output of a 'call' from its fully qualified name (e.g. my_workflow.my_call).",
    nickname = "call-stdout",
    httpMethod = "GET",
    produces = "text/plain"
    //,response = classOf[CallOutputs]
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "version", required = true, dataType = "string", paramType = "path", value = "API Version", allowableValues = CromwellApiService.VersionAllowableValues),
    new ApiImplicitParam(name = "workflowId", required = true, dataType = "string", paramType = "path", value = "Workflow ID"),
    new ApiImplicitParam(name = "callFqn", required = true, dataType = "string", paramType = "path", value = "Call fully qualified name")
  ))
  @ApiResponses(Array(
    // See the note above regarding Swagger brokenness for why there is no response classOf here.
    new ApiResponse(code = 200, message = "Successful Request"),
    new ApiResponse(code = 404, message = "Workflow ID Not Found"),
    new ApiResponse(code = 404, message = "Call Fully Qualified Name Not Found"),
    new ApiResponse(code = 400, message = "Malformed Workflow ID"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def callStdoutStderrRoute =
    path("workflows" / Segment / Segment / "logs" / Segment) { (version, workflowId, callFqn) =>
      Try(UUID.fromString(workflowId)) match {
        case Success(w) =>
          // This currently does not attempt to parse the call name for conformation to any pattern.
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.CallStdoutStderr(w, callFqn))
        case Failure(ex) =>
          complete(StatusCodes.BadRequest, s"Invalid workflow ID: '$workflowId'.")
      }
    }
}
