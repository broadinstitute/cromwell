package cromwell.docker.registryv2.flows.ecr

import cats.effect.IO
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.{DockerImageIdentifier, DockerInfoActor, DockerRegistryConfig}
import org.http4s.client.Client
import software.amazon.awssdk.services.ecr.EcrClient

import scala.collection.JavaConverters._

class EcrRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {

  override def getToken(dockerInfoContext: DockerInfoActor.DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {
    IO(EcrClient.create().getAuthorizationToken().authorizationData().asScala.headOption.map(_.authorizationToken()))
  }

  override def accepts(dockerImageIdentifier: DockerImageIdentifier) = {
    dockerImageIdentifier.hostAsString.contains(".dkr.ecr.")
  }
}
