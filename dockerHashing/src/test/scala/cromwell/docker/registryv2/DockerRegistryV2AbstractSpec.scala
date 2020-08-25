package cromwell.docker.registryv2

import cats.effect.{IO, Resource}
import cromwell.docker.DockerInfoActor.{DockerInfoContext, DockerInfoFailedResponse}
import cromwell.docker.{DockerImageIdentifier, DockerInfoActor, DockerInfoRequest, DockerRegistryConfig}
import org.http4s.client.Client
import org.http4s.headers.`Content-Type`
import org.http4s.{Header, Headers, MediaType, Request, Response}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class DockerRegistryV2AbstractSpec extends AnyFlatSpec with Matchers {
  behavior of "DockerRegistryV2Abstract"
  
  it should "handle gracefully if a response cannot be parsed" in {
    val registry = new DockerRegistryV2Abstract(DockerRegistryConfig.default) {
      override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier) = "N/A"
      override protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier) = "N/A"
      override protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoActor.DockerInfoContext) = List.empty
    }
    
    val mediaType = MediaType.parse(DockerRegistryV2Abstract.ManifestV2MediaType).right.get
    val contentType: Header = `Content-Type`(mediaType)

    val mockClient = Client({ _: Request[IO] =>
      // This response will have an empty body, so we need to be explicit about the typing:
      Resource.pure[IO, Response[IO]](Response(headers = Headers.of(contentType))) : Resource[IO, Response[IO]]
    })
    
    val dockerImageIdentifier = DockerImageIdentifier.fromString("ubuntu").get
    val dockerInfoRequest = DockerInfoRequest(dockerImageIdentifier)
    val context = DockerInfoContext(dockerInfoRequest, null)
    val result = registry.run(context)(mockClient).unsafeRunSync()
    result.asInstanceOf[(DockerInfoFailedResponse, DockerInfoContext)]._1.reason shouldBe "Failed to get docker hash for ubuntu:latest Malformed message body: Invalid JSON: empty body"
  }
}
