package cromwell.docker.local

import cromwell.core.Tags.IntegrationTest
import cromwell.docker.DockerInfoActor.{DockerInfoNotFound, DockerInfoSuccessResponse, DockerInformation}
import cromwell.docker.{DockerHashResult, DockerRegistry, DockerRegistrySpec}
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._

class DockerCliSpec extends DockerRegistrySpec("DockerCliFlowSpec") with FlatSpecLike with Matchers {
  behavior of "DockerCliFlow"

  override protected def registryFlows: Seq[DockerRegistry] = Seq(new DockerCliFlow)

  it should "retrieve a public docker hash" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("ubuntu:latest")

    expectMsgPF(30.seconds) {
      case DockerInfoSuccessResponse(DockerInformation(DockerHashResult(alg, hash), _), _) =>
        alg shouldBe "sha256"
        hash should not be empty
    }
  }

  it should "retrieve a public docker hash on gcr" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("gcr.io/google-containers/alpine-with-bash:1.0")

    expectMsgPF(30.seconds) {
      case DockerInfoSuccessResponse(DockerInformation(DockerHashResult(alg, hash), _), _) =>
        alg shouldBe "sha256"
        hash should not be empty
    }
  }

  it should "send image not found message back if the image does not exist" taggedAs IntegrationTest in {
    val notFound = makeRequest("ubuntu:nonexistingtag")
    dockerActor ! notFound

    expectMsgClass(30.seconds, classOf[DockerInfoNotFound])
  }

  it should "send image not found if the user doesn't have permission on this image" taggedAs IntegrationTest in {
    val unauthorized = makeRequest("tjeandet/sinatra:v1")
    dockerActor ! unauthorized

    expectMsgClass(30.seconds, classOf[DockerInfoNotFound])
  }

  it should "send image not found for an unrecognized host" taggedAs IntegrationTest in {
    val unauthorized = makeRequest("unknown.io/image:v1")
    dockerActor ! unauthorized

    expectMsgClass(30.seconds, classOf[DockerInfoNotFound])
  }
}
