package cromwell.webservice

import java.util.UUID

import akka.actor.{Actor, ActorRef, ActorRefFactory, Props}
import com.gettyimages.spray.swagger.SwaggerHttpService
import com.wordnik.swagger.annotations._
import com.wordnik.swagger.model.ApiInfo
import spray.http.StatusCodes
import spray.json._
import spray.httpx.SprayJsonSupport._
import spray.routing.Directive.pimpApply
import spray.routing._

import scala.reflect.runtime.universe._
import scala.util.{Failure, Success, Try}


object CromwellApiServiceActor {
  def props(workflowManagerActorRef : ActorRef, swaggerService: SwaggerService): Props = {
    Props(new CromwellApiServiceActor(workflowManagerActorRef, swaggerService))
  }
}

class SwaggerService(override val apiVersion: String,
                     override val baseUrl: String,
                     override val docsPath: String,
                     override val swaggerVersion: String,
                     override val apiTypes: Seq[Type],
                     override val apiInfo: Option[ApiInfo])
                    (implicit val actorRefFactory: ActorRefFactory)
  extends SwaggerHttpService

class CromwellApiServiceActor(val workflowManager : ActorRef, swaggerService: SwaggerService) extends Actor with CromwellApiService {
  implicit def executionContext = actorRefFactory.dispatcher
  def actorRefFactory = context

  def possibleRoutes = workflowRoutes ~ swaggerService.routes ~
    get {
      pathSingleSlash {
        getFromResource("swagger/index.html")
      } ~ getFromResourceDirectory("swagger/") ~ getFromResourceDirectory("META-INF/resources/webjars/swagger-ui/2.0.24/")
    }

  def receive = runRoute(possibleRoutes)

}

@Api(value = "/workflow", description = "Workflow", produces = "application/json", position = 1)
trait CromwellApiService extends HttpService with PerRequestCreator  {
  val workflowManager: ActorRef

  val workflowRoutes = queryRoute ~ outputsRoute ~ submitRoute

  @ApiOperation(
    value = "Query for workflow status based on workflow id.",
    nickname = "status",
    httpMethod = "GET",
    produces = "application/json"
    //, response = classOf[WorkflowStatusResponse]
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "id", required = true, dataType = "string", paramType = "path", value = "workflow identifier")
  ))
  @ApiResponses(Array(
    new ApiResponse(code = 200, message = "Successful Request"),
    new ApiResponse(code = 404, message = "Workflow ID Not Found"),
    new ApiResponse(code = 400, message = "Malformed Workflow ID"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def queryRoute =
    path("workflow" / Segment / "status") { id =>
      Try(UUID.fromString(id)) match {
        case Success(workflowId) =>
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowStatus(workflowId))
        case Failure(ex) =>
          complete(StatusCodes.BadRequest)
      }
    }

  @ApiOperation(
    value = "Submit a new workflow for execution",
    nickname = "submit",
    httpMethod = "POST",
    produces = "application/json"
//  ,response = classOf[WorkflowSubmitResponse]
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "wdlSource", required = true, dataType = "string", paramType = "form"),
    new ApiImplicitParam(name = "workflowInputs", required = true, dataType = "string", paramType = "form")
  ))
  @ApiResponses(Array(
    new ApiResponse(code = 201, message = "Successful Request"),
    new ApiResponse(code = 400, message = "Malformed Input"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def submitRoute =
    path("workflows") {
      post {
        formFields("wdlSource", "workflowInputs") { (wdlSource, workflowInputs) =>

          val tryMap = Try(workflowInputs.parseJson)
          tryMap match {
            case Success(JsObject(json)) =>
              requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.SubmitWorkflow(wdlSource, workflowInputs, json))
            case Success(o) =>
              complete(StatusCodes.BadRequest, "Expecting JSON object as workflow inputs")
            case Failure(ex) =>
              complete(StatusCodes.BadRequest, "workflowInput JSON was malformed")
          }
        }
      }
    }

  @ApiOperation(
    value = "Query for workflow outputs based on workflow id.",
    nickname = "outputs",
    httpMethod = "GET",
    produces = "application/json"
    //,response = classOf[WorkflowOutputs]
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "id", required = true, dataType = "string", paramType = "path", value = "workflow")
  ))
  @ApiResponses(Array(
    new ApiResponse(code = 200, message = "Successful Request"),
    new ApiResponse(code = 404, message = "Workflow ID Not Found"),
    new ApiResponse(code = 400, message = "Malformed Workflow ID"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def outputsRoute =
    path("workflow" / Segment / "outputs") { id =>
      Try(UUID.fromString(id)) match {
        case Success(workflowId) =>
          requestContext => perRequest(requestContext, CromwellApiHandler.props(workflowManager), CromwellApiHandler.WorkflowOutputs(workflowId))
        case Failure(ex) =>
          complete(StatusCodes.BadRequest)
    }
  }
}
