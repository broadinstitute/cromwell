package cromwell.docker.registryv2.flows.aws

import cats.effect.IO
import cromwell.docker.{DockerImageIdentifier, DockerInfoActor, DockerRegistryConfig}
import org.http4s.client.Client
import software.amazon.awssdk.services.ecrpublic.EcrPublicClient
import software.amazon.awssdk.services.ecrpublic.model.GetAuthorizationTokenRequest

import scala.concurrent.Future


class AmazonEcrPublic(override val config: DockerRegistryConfig, ecrClient: EcrPublicClient = EcrPublicClient.create()) extends AmazonEcrAbstract(config) {
  /**
    * public.ecr.aws
    */
  override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String = "public.ecr.aws"

  /**
    * Returns true if this flow is able to process this docker image,
    * false otherwise
    */
  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean = dockerImageIdentifier.hostAsString.contains("public.ecr.aws")


  override protected def getToken(dockerInfoContext: DockerInfoActor.DockerInfoContext)(implicit client: Client[IO]): IO[Option[String]] = {

    val eventualMaybeToken: Future[Option[String]] = Future(
      Option(ecrClient
        .getAuthorizationToken(GetAuthorizationTokenRequest.builder().build())
        .authorizationData.authorizationToken()
      )
    )

    IO.fromFuture(IO(eventualMaybeToken))
  }
}
