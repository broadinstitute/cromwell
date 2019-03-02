package cromiam.webservice

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import akka.http.scaladsl.model.StatusCodes._
import akka.http.scaladsl.model.{HttpResponse, StatusCodes}
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import cromiam.cromwell.CromwellClient
import cromiam.server.status.StatusService
import cromiam.webservice.EngineRouteSupport._
import cromwell.api.model._
import org.broadinstitute.dsde.workbench.util.health.StatusJsonSupport._

import scala.concurrent.ExecutionContextExecutor


trait EngineRouteSupport extends RequestSupport with SprayJsonSupport {
  val statusService: StatusService
  val cromwellClient: CromwellClient

  implicit def executor: ExecutionContextExecutor

  val engineRoutes: Route = statsRoute ~ versionRoute ~ engineStatusRoute

  def versionRoute: Route = path("engine" / Segment / "version") { _ =>
    get {
      extractStrictRequest { req =>
        complete { cromwellClient.forwardToCromwell(req).asHttpResponse }
      }
    }
  }

  def engineStatusRoute: Route = path("engine" / Segment / "status") { _ =>
    get {
      complete(statusService.status().map { statusResponse =>
        val httpStatus = if (statusResponse.ok) StatusCodes.OK else StatusCodes.InternalServerError
        (httpStatus, statusResponse)
      })
    }
  }

  def statsRoute: Route = path("engine" / Segment / "stats") { _ => complete(CromIamStatsForbidden) }

}

object EngineRouteSupport {
  private[webservice] val CromIamStatsForbidden = HttpResponse(status = Forbidden, entity = "CromIAM does not allow access to the /stats endpoint")
}
