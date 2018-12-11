package cromwell.docker

import cromwell.docker.DockerInfoActor.DockerHashUnknownRegistry
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._

class DockerEmptyFlowSpec extends DockerRegistrySpec("DockerEmptyFlowSpec") with FlatSpecLike with Matchers {
  behavior of "An empty docker flow"

  override protected def registryFlows: Seq[DockerRegistry] = Seq()

  it should "send an unrecognized host message back for a public docker hash" in {
    dockerActor ! makeRequest("ubuntu:latest")

    expectMsgClass(5.seconds, classOf[DockerHashUnknownRegistry])
  }

  it should "send an unrecognized host message back for a public docker on gcr" in {
    dockerActor ! makeRequest("gcr.io/google-containers/alpine-with-bash:1.0")

    expectMsgClass(5.seconds, classOf[DockerHashUnknownRegistry])
  }

  it should "send an unrecognized host message back if the image does not exist" in {
    val notFound = makeRequest("ubuntu:nonexistingtag")
    dockerActor ! notFound

    expectMsgClass(5.seconds, classOf[DockerHashUnknownRegistry])
  }

  it should "send an unrecognized host message back if the user doesn't have permission on this image" in {
    val unauthorized = makeRequest("tjeandet/sinatra:v1")
    dockerActor ! unauthorized

    expectMsgClass(5.seconds, classOf[DockerHashUnknownRegistry])
  }

  it should "send an unrecognized host message back for an unrecognized host" in {
    val unauthorized = makeRequest("unknown.io/image:v1")
    dockerActor ! unauthorized

    expectMsgClass(5.seconds, classOf[DockerHashUnknownRegistry])
  }
}
