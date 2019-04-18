package cromwell.docker.registryv2.flows.quay

import cats.effect.IO
import cromwell.docker.DockerInfoActor.DockerInfoContext
import cromwell.docker.DockerRegistryConfig
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import org.http4s.client.Client

class QuayRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {
  override protected def getToken(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {
    IO.pure(None)
  }
}
