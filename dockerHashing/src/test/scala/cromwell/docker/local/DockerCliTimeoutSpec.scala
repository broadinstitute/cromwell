package cromwell.docker.local

import cromwell.core.Tags.IntegrationTest
import cromwell.docker.DockerInfoActor.DockerInfoFailedResponse
import cromwell.docker.{DockerRegistry, DockerRegistrySpec}
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.TimeoutException
import scala.concurrent.duration._

class DockerCliTimeoutSpec extends DockerRegistrySpec("DockerCliTimeoutFlowSpec") with FlatSpecLike with Matchers {
  behavior of "A DockerCliFlow that times out"

  override protected def registryFlows: Seq[DockerRegistry] = Seq(new DockerCliFlow {
    override lazy val firstLookupTimeout = 0.seconds
  })

  it should "timeout retrieving a public docker hash" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("ubuntu:latest")

    expectMsgPF(5.seconds) {
      case DockerInfoFailedResponse(exception: TimeoutException, _) =>
        exception.getMessage should be(
          """|Timeout while looking up hash of ubuntu:latest.
             |Ensure that docker is running correctly.
             |""".stripMargin)
    }
  }

  it should "timeout retrieving a public docker hash on gcr" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("gcr.io/google-containers/alpine-with-bash:1.0")

    expectMsgPF(5.seconds) {
      case DockerInfoFailedResponse(exception: TimeoutException, _) =>
        exception.getMessage should be(
          """|Timeout while looking up hash of gcr.io/google-containers/alpine-with-bash:1.0.
             |Ensure that docker is running correctly.
             |""".stripMargin)
    }

  }

  it should "timeout retrieving an image that does not exist" taggedAs IntegrationTest in {
    val notFound = makeRequest("ubuntu:nonexistingtag")
    dockerActor ! notFound

    expectMsgPF(5.seconds) {
      case DockerInfoFailedResponse(exception: TimeoutException, _) =>
        exception.getMessage should be(
          """|Timeout while looking up hash of ubuntu:nonexistingtag.
             |Ensure that docker is running correctly.
             |""".stripMargin)
    }
  }

  it should "timeout retrieving an image if the user doesn't have permission" taggedAs IntegrationTest in {
    val unauthorized = makeRequest("tjeandet/sinatra:v1")
    dockerActor ! unauthorized

    expectMsgPF(5.seconds) {
      case DockerInfoFailedResponse(exception: TimeoutException, _) =>
        exception.getMessage should be(
          """|Timeout while looking up hash of tjeandet/sinatra:v1.
             |Ensure that docker is running correctly.
             |""".stripMargin)
    }
  }

  it should "timeout retrieving an image for an unrecognized host" taggedAs IntegrationTest in {
    val unauthorized = makeRequest("unknown.io/image:v1")
    dockerActor ! unauthorized

    expectMsgPF(5.seconds) {
      case DockerInfoFailedResponse(exception: TimeoutException, _) =>
        exception.getMessage should be(
          """|Timeout while looking up hash of unknown.io/library/image:v1.
             |Ensure that docker is running correctly.
             |""".stripMargin)
    }
  }
}
