package cromwell.webservice

import java.util.UUID

import akka.actor.{ActorRef, Actor, ActorRefFactory, Props}
import com.gettyimages.spray.swagger.SwaggerHttpService
import com.wordnik.swagger.annotations._
import com.wordnik.swagger.model.ApiInfo
import cromwell.engine.WorkflowId
import spray.http.StatusCodes
import spray.routing.Directive.pimpApply
import spray.routing._
import scala.util.{Success, Failure}

import scala.reflect.runtime.universe._
import scala.util.Try

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

class CromwellApiServiceActor(val workflowManagerActorRef : ActorRef, swaggerService: SwaggerService) extends Actor with CromwellApiService {
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

@Api(value = "/workflows", description = "Workflow", produces = "application/json", position = 1)
trait CromwellApiService extends HttpService with PerRequestCreator  {
  val workflowManagerActorRef : ActorRef

  val workflowRoutes = queryRoute

  @ApiOperation(value = "Query for workflow based on workflow id.",
    nickname = "query",
    httpMethod = "GET",
    produces = "application/json",
    response = classOf[WorkflowStatusResponse]
  )
  @ApiImplicitParams(Array(
    new ApiImplicitParam(name = "id", required = true, dataType = "string", paramType = "path", value = "workflow")
  ))
  @ApiResponses(Array(
    new ApiResponse(code = 200, message = "Successful Request"),
    new ApiResponse(code = 404, message = "Workflow ID Not Found"),
    new ApiResponse(code = 500, message = "Internal Error")
  ))
  def queryRoute =
    path("workflows" / Segment) { id =>
      Try(UUID.fromString(id)) match {
        case Success(workflowId) =>
          requestContext =>
            perRequest(requestContext, CromwellApiHandler.props(workflowManagerActorRef), CromwellApiHandler.WorkflowStatus(workflowId))
        case Failure(ex) =>
          complete(StatusCodes.BadRequest)
      }
    }
}
