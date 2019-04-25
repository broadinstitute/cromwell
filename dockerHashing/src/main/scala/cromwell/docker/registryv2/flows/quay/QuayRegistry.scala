package cromwell.docker.registryv2.flows.quay

import cats.effect.IO
import cromwell.docker.DockerInfoActor.DockerInfoContext
import cromwell.docker.{DockerImageIdentifier, DockerRegistryConfig}
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import org.http4s.client.Client

class QuayRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {
  override protected def getToken(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {
    IO.pure(None)
  }

  /**
    * Returns true if this flow is able to process this docker image,
    * false otherwise
    */
  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean = {
    dockerImageIdentifier.hostAsString.startsWith("quay.io")
  }
}
