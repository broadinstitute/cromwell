package cromwell.docker.registryv2.flows.quay

import cats.effect.IO
import cromwell.docker.DockerInfoActor.DockerInfoContext
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.{DockerImageIdentifier, DockerRegistryConfig}
import org.http4s.Header
import org.http4s.client.Client

class QuayRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {
  override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String = "quay.io"
  // Not used for now because we bypass the token part as quay doesn't require one for public images 
  override protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier): String = "quay.io"
  // Not used for now, same reason as above
  override protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext): List[Header] = List.empty

  override protected def getToken(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {
    IO.pure(None)
  }
}
