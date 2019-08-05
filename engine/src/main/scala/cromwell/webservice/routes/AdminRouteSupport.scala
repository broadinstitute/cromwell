package cromwell.webservice.routes

import akka.actor.{ActorRef, ActorRefFactory}
import akka.http.scaladsl.server.Directives._
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.webservice.WebServiceUtils

import scala.concurrent.ExecutionContext

trait AdminRouteSupport extends WebServiceUtils {

  implicit def actorRefFactory: ActorRefFactory
  implicit val ec: ExecutionContext
  implicit val materializer: ActorMaterializer
  implicit val timeout: Timeout

  val serviceRegistryActor: ActorRef

  val adminRoutes = concat(
    path("admin" / Segment / "listSubmission") { _ => complete("listSubmissions") },
    path("admin" / Segment / "pauseSubmission") { _ => complete("pauseSubmissions") }
  )
}
