package cromwell.docker.registryv2

import cats.effect.{IO, Resource}
import cromwell.docker.DockerInfoActor.{DockerInfoContext, DockerInfoFailedResponse}
import cromwell.docker.{DockerImageIdentifier, DockerInfoActor, DockerInfoRequest, DockerRegistryConfig}
import org.http4s.client.Client
import org.http4s.headers.`Content-Type`
import org.http4s.{Header, Headers, MediaType, Response}
import org.scalatest.{FlatSpec, Matchers}

class DockerRegistryV2AbstractSpec extends FlatSpec with Matchers {
  behavior of "DockerRegistryV2Abstract"
  
  it should "handle gracefully if a response cannot be parsed" in {
    val registry = new DockerRegistryV2Abstract(DockerRegistryConfig.default) {
      override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier) = "N/A"
      override protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier) = "N/A"
      override protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoActor.DockerInfoContext) = List.empty
    }
    
    val mediaType = MediaType.parse(DockerRegistryV2Abstract.ManifestV2MediaType).right.get
    val contentType: Header = `Content-Type`(mediaType)
    
    val mockClient = Client[IO]({ _ =>
      // This response will have an empty body
      Resource.pure(
        Response(headers = Headers(contentType))
      )
    })
    
    val dockerImageIdentifier = DockerImageIdentifier.fromString("ubuntu").get
    val dockerInfoRequest = DockerInfoRequest(dockerImageIdentifier)
    val context = DockerInfoContext(dockerInfoRequest, null)
    val result = registry.run(context)(mockClient).unsafeRunSync()
    result.asInstanceOf[(DockerInfoFailedResponse, DockerInfoContext)]._1.reason shouldBe "Failed to get docker hash for ubuntu:latest Malformed message body: Invalid JSON: empty body"
  }
}
