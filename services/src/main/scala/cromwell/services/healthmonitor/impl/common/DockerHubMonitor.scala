package cromwell.services.healthmonitor.impl.common

import akka.actor.ActorSystem
import akka.pattern.ask
import akka.util.Timeout
import cats.effect.{ContextShift, IO}
import cromwell.docker.DockerInfoActor.{DockerInfoFailedResponse, DockerInfoResponse, DockerInfoSuccessResponse}
import cromwell.docker.registryv2.flows.dockerhub.DockerHubRegistry
import cromwell.docker.{DockerImageIdentifierWithoutHash, DockerInfoActor, DockerInfoRequest, DockerRegistryConfig}
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}

/**
  * A mixin which provides a dockerhub connectivity monitor. Ideally this would be an Object and not a Trait but the
  * need for an ExecutionContext made that difficult with the () => Future[Subsystem] signature of the checks.
  */
trait DockerHubMonitor {
  import DockerHubMonitor.UbuntuLatestHashRequest
  implicit val ec: ExecutionContext
  implicit lazy val cs: ContextShift[IO] = IO.contextShift(ec)
  implicit val system: ActorSystem
  implicit val timeout = Timeout(5.seconds)

  lazy val dockerHubFlow = List(new DockerHubRegistry(DockerRegistryConfig.default))
  lazy val dockerHashActor = system.actorOf(DockerInfoActor.props(dockerHubFlow, 500, 0.minutes, 0), "HealthMonitorDockerHashActor")

  lazy val DockerHub = MonitoredSubsystem("DockerHub", checkDockerhub _)

  /**
    * Demonstrates connectivity to Docker Hub by periodically pulling the hash of an image. If the hash is not returned
    * we assume that there is no connection to Docker Hub.
    */
  private def checkDockerhub(): Future[SubsystemStatus] = {
    dockerHashActor.ask(UbuntuLatestHashRequest).mapTo[DockerInfoResponse] map {
      case _: DockerInfoSuccessResponse => OkStatus
      case f: DockerInfoFailedResponse => throw f.failure
      case huh => throw new RuntimeException("Encountered unexpected error when trying to contact DockerHub: " + huh.getClass.getCanonicalName)
    }
  }
}

object DockerHubMonitor {
  val UbuntuLatestHashRequest = DockerInfoRequest(DockerImageIdentifierWithoutHash(None, None, "ubuntu", "latest"))
}
