package cromwell.docker.registryv2.flows.aws

import cats.effect.IO
import cromwell.docker.{DockerImageIdentifier, DockerInfoActor, DockerRegistryConfig}
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.registryv2.flows.aws.AwsElasticContainerRegistry.{isEcr, isPublicEcr}
import org.http4s.{AuthScheme, Header}
import org.http4s.client.Client
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.ecr.EcrClient
import software.amazon.awssdk.services.ecrpublic.EcrPublicClient
import software.amazon.awssdk.services.ecrpublic.model.GetAuthorizationTokenRequest

import scala.compat.java8.OptionConverters.RichOptionalGeneric

class AwsElasticContainerRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {

  private lazy val ecrClient = EcrClient.create()
  private lazy val ecrPublicClient = EcrPublicClient.builder().region(Region.US_EAST_1).build()

  override def getAuthorizationScheme(dockerImageIdentifier: DockerImageIdentifier): AuthScheme =
    if (isPublicEcr(dockerImageIdentifier.hostAsString)) AuthScheme.Bearer else AuthScheme.Basic

  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean =
    isEcr(dockerImageIdentifier.hostAsString)

  /**
   * (e.g registry-1.docker.io)
   */
  override def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String =
    dockerImageIdentifier.host.getOrElse("")

  /**
   * (e.g auth.docker.io)
   */
  override def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier): String =
    dockerImageIdentifier.host.getOrElse("")

  override def getToken(
    dockerInfoContext: DockerInfoActor.DockerInfoContext
  )(implicit client: Client[IO]): IO[Option[String]] =
    if (isPublicEcr(dockerInfoContext.dockerImageID.hostAsString)) getPublicEcrToken
    else getPrivateEcrToken

  /**
   * Builds the list of headers for the token request
   */
  override def buildTokenRequestHeaders(dockerInfoContext: DockerInfoActor.DockerInfoContext): List[Header] =
    List.empty

  private def getPublicEcrToken: IO[Option[String]] =
    IO(
      Option(
        ecrPublicClient
          .getAuthorizationToken(GetAuthorizationTokenRequest.builder().build())
          .authorizationData()
          .authorizationToken()
      )
    )

  private def getPrivateEcrToken: IO[Option[String]] =
    IO(
      ecrClient
        .getAuthorizationToken()
        .authorizationData()
        .stream()
        .findFirst()
        .asScala
        .map(_.authorizationToken())
    )

}

object AwsElasticContainerRegistry {
  def isEcr(host: String): Boolean = isPublicEcr(host) || isPrivateEcr(host)
  def isPublicEcr(host: String): Boolean = host.contains("public.ecr.aws")
  def isPrivateEcr(host: String): Boolean = host.contains("amazonaws.com")
}
