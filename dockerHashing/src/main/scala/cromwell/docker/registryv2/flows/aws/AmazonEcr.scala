package cromwell.docker.registryv2.flows.aws

import cats.effect.IO
import cromwell.docker.{DockerImageIdentifier, DockerInfoActor, DockerRegistryConfig}
import org.http4s.AuthScheme
import org.http4s.client.Client
import software.amazon.awssdk.services.ecr.EcrClient

import scala.compat.java8.OptionConverters._
import scala.concurrent.Future

class AmazonEcr(override val config: DockerRegistryConfig, ecrClient: EcrClient = EcrClient.create()) extends AmazonEcrAbstract(config) {

  override protected val authorizationScheme: AuthScheme = AuthScheme.Basic

  /**
    * e.g 123456789012.dkr.ecr.us-east-1.amazonaws.com
    */
  override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String = {
    var hostname = dockerImageIdentifier.hostAsString
    if (hostname.lastIndexOf("/").equals(hostname.length -1)) {
      hostname = hostname.substring(0, hostname.length -1)
    }
    hostname
  }
  /**
    * Returns true if this flow is able to process this docker image,
    * false otherwise
    */
  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean = dockerImageIdentifier.hostAsString.contains("amazonaws.com")

  override protected def getToken(dockerInfoContext: DockerInfoActor.DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {
    val eventualMaybeToken = Future(ecrClient.getAuthorizationToken
      .authorizationData()
      .stream()
      .findFirst()
      .asScala
      .map(_.authorizationToken()))

    IO.fromFuture(IO(eventualMaybeToken))
  }
}
