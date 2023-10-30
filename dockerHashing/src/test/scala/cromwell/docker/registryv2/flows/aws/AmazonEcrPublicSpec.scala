package cromwell.docker.registryv2.flows.aws
import cats.effect.{IO, Resource}
import cromwell.core.TestKitSuite
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.{DockerImageIdentifier, DockerInfoActor, DockerInfoRequest, DockerRegistryConfig}
import org.http4s.{Header, Headers, MediaType, Request, Response}
import org.http4s.client.Client
import org.http4s.headers.`Content-Type`
import org.mockito.ArgumentMatchers.any
import org.mockito.Mockito._
import org.scalatest.{BeforeAndAfter, PrivateMethodTester}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import software.amazon.awssdk.services.ecrpublic.model.{AuthorizationData, GetAuthorizationTokenRequest, GetAuthorizationTokenResponse}
import software.amazon.awssdk.services.ecrpublic.EcrPublicClient

class AmazonEcrPublicSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with BeforeAndAfter with PrivateMethodTester {
  behavior of "AmazonEcrPublic"

  val goodUri = "public.ecr.aws/amazonlinux/amazonlinux:latest"
  val otherUri = "ubuntu:latest"


  val mediaType: MediaType = MediaType.parse(DockerRegistryV2Abstract.DockerManifestV2MediaType).getOrElse(fail("Cant parse mediatype"))
  val contentType: Header = `Content-Type`(mediaType)
  val mockEcrClient: EcrPublicClient = mock(classOf[EcrPublicClient])
  implicit val mockIOClient: Client[IO] = Client({ _: Request[IO] =>
    // This response will have an empty body, so we need to be explicit about the typing:
    Resource.pure[IO, Response[IO]](Response(headers = Headers.of(contentType))) : Resource[IO, Response[IO]]
  })

  val registry = new AmazonEcrPublic(DockerRegistryConfig.default, mockEcrClient)

  it should "Accept good URI" in {
    val dockerImageIdentifier = DockerImageIdentifier.fromString(goodUri).get
    registry.accepts(dockerImageIdentifier) shouldEqual true
  }

  it should "NOT accept other URI" in {
    val dockerImageIdentifier = DockerImageIdentifier.fromString(otherUri).get
    registry.accepts(dockerImageIdentifier) shouldEqual false
  }

  it should "have public.ecr.aws as registryHostName" in {
    val registryHostNameMethod = PrivateMethod[String](Symbol("registryHostName"))
    registry invokePrivate registryHostNameMethod(DockerImageIdentifier.fromString(goodUri).get) shouldEqual "public.ecr.aws"
  }

  it should "return expected auth token" in {
    val token = "auth-token"
    val imageId = DockerImageIdentifier.fromString(goodUri).get
    val dockerInfoRequest = DockerInfoRequest(imageId)
    val context = DockerInfoActor.DockerInfoContext(request = dockerInfoRequest, replyTo =  emptyActor)

    when(mockEcrClient.getAuthorizationToken(any[GetAuthorizationTokenRequest]()))
      .thenReturn(GetAuthorizationTokenResponse
        .builder()
        .authorizationData(AuthorizationData
          .builder()
          .authorizationToken(token)
          .build())
        .build)

    val getTokenMethod = PrivateMethod[IO[Option[String]]](Symbol("getToken"))
    registry invokePrivate getTokenMethod(context, mockIOClient) ensuring(io => io.unsafeRunSync().get == token)
  }
}
