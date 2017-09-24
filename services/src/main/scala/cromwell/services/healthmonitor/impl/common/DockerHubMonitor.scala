package cromwell.services.healthmonitor.impl.common

import akka.actor.ActorSystem
import akka.pattern.ask
import akka.http.scaladsl.Http
import akka.stream.ActorMaterializer
import akka.util.Timeout
import cromwell.docker.DockerHashActor.{DockerHashContext, DockerHashFailedResponse, DockerHashResponse, DockerHashSuccessResponse}
import cromwell.docker.registryv2.flows.HttpFlowWithRetry.ContextWithRequest
import cromwell.docker.registryv2.flows.dockerhub.DockerHubFlow
import cromwell.docker.{DockerHashActor, DockerHashRequest, DockerImageIdentifierWithoutHash}
import cromwell.services.healthmonitor.HealthMonitorServiceActor._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

/**
  * A mixin which provides a dockerhub connectivity monitor. Ideally this would be an Object and not a Trait but the
  * need for an ExecutionContext made that difficult with the () => Future[Subsystem] signature of the checks.
  */
trait DockerHubMonitor {
  import DockerHubMonitor.UbuntuLatestHashRequest
  implicit val ec: ExecutionContext
  implicit val system: ActorSystem
  implicit lazy val materializer = ActorMaterializer()
  implicit val timeout = Timeout(5.seconds)

  lazy val dockerHttpPool = Http().superPool[ContextWithRequest[DockerHashContext]]()
  lazy val dockerHubFlow = List(new DockerHubFlow(dockerHttpPool)(ec, materializer, system.scheduler))
  lazy val dockerHashActor = system.actorOf(DockerHashActor.props(dockerHubFlow, 500, 0.minutes, 0)(materializer), "DockerHashActor")

  lazy val DockerHub = MonitoredSubsystem("DockerHub", checkDockerhub _)

  /**
    * Demonstrates connectivity to DockerHub by periodically pulling a small image. Remove the image afterwards in
    * order to ensure that subsequent attempts aren't showing a false positive
    */
  private def checkDockerhub(): Future[SubsystemStatus] = {
    dockerHashActor.ask(UbuntuLatestHashRequest).mapTo[DockerHashResponse] map {
      case _: DockerHashSuccessResponse => OkStatus
      case f: DockerHashFailedResponse => throw f.failure
      case huh => throw new RuntimeException("Encountered unexpected error when trying to contact DockerHub: " + huh.getClass.getCanonicalName)
    }
  }
}

object DockerHubMonitor {
  val UbuntuLatestHashRequest = DockerHashRequest(DockerImageIdentifierWithoutHash(None, None, "ubuntu", "latest"))
}
