package cromwell.docker.registryv2.flows.aws

import cats.effect.{IO, Resource}
import cromwell.core.TestKitSuite
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.{DockerImageIdentifier, DockerInfoActor, DockerInfoRequest, DockerRegistryConfig}
import org.http4s.{AuthScheme, Header, Headers, MediaType, Request, Response}
import org.http4s.client.Client
import org.http4s.headers.`Content-Type`
import org.mockito.Mockito._
import org.scalatest.{BeforeAndAfter, PrivateMethodTester}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatestplus.mockito.MockitoSugar
import software.amazon.awssdk.services.ecr.EcrClient
import software.amazon.awssdk.services.ecr.model.{AuthorizationData, GetAuthorizationTokenResponse}

class AmazonEcrSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with MockitoSugar with BeforeAndAfter with PrivateMethodTester{
  behavior of "AmazonEcr"

  val goodUri = "123456789012.dkr.ecr.us-east-1.amazonaws.com/amazonlinux/amazonlinux:latest"
  val otherUri = "ubuntu:latest"

  val mediaType: MediaType = MediaType.parse(DockerRegistryV2Abstract.ManifestV2MediaType).right.get
  val contentType: Header = `Content-Type`(mediaType)
  val mockEcrClient: EcrClient = mock[EcrClient]
  implicit val mockIOClient: Client[IO] = Client({ _: Request[IO] =>
    // This response will have an empty body, so we need to be explicit about the typing:
    Resource.pure[IO, Response[IO]](Response(headers = Headers.of(contentType))) : Resource[IO, Response[IO]]
  })

  val registry = new AmazonEcr(DockerRegistryConfig.default, mockEcrClient)

  it should "accept good URI" in {
    val dockerImageIdentifier = DockerImageIdentifier.fromString(goodUri).get
    registry.accepts(dockerImageIdentifier) shouldEqual true
  }

  it should "NOT accept other URI" in {
    val dockerImageIdentifier = DockerImageIdentifier.fromString(otherUri).get
    registry.accepts(dockerImageIdentifier) shouldEqual false
  }

  it should "use Basic Auth Scheme" in {
    val authSchemeMethod = PrivateMethod[AuthScheme]('authorizationScheme)
    registry invokePrivate authSchemeMethod() shouldEqual AuthScheme.Basic
  }

  it should "return 123456789012.dkr.ecr.us-east-1.amazonaws.com as registryHostName" in {
    val registryHostNameMethod = PrivateMethod[String]('registryHostName)
    registry invokePrivate registryHostNameMethod(DockerImageIdentifier.fromString(goodUri).get) shouldEqual "123456789012.dkr.ecr.us-east-1.amazonaws.com"
  }

  it should "return expected auth token" in {
    val token = "auth-token"
    val imageId = DockerImageIdentifier.fromString(goodUri).get
    val dockerInfoRequest = DockerInfoRequest(imageId)
    val context = DockerInfoActor.DockerInfoContext(request = dockerInfoRequest, replyTo =  emptyActor)

    when(mockEcrClient.getAuthorizationToken)
      .thenReturn(GetAuthorizationTokenResponse
        .builder()
        .authorizationData(AuthorizationData
          .builder()
          .authorizationToken(token)
          .build())
        .build)

    val getTokenMethod = PrivateMethod[IO[Option[String]]]('getToken)
    registry invokePrivate getTokenMethod(context, mockIOClient) ensuring(io => io.unsafeRunSync().get == token)
  }
}
