package cromwell.docker

import cromwell.core.Tags.IntegrationTest
import cromwell.docker.DockerInfoActor._
import cromwell.docker.registryv2.flows.dockerhub.DockerHubRegistry
import cromwell.docker.registryv2.flows.gcr.GcrRegistry
import cromwell.docker.registryv2.flows.quay.QuayRegistry
import org.scalatest.{BeforeAndAfterAll, FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class DockerInfoActorSpec extends DockerRegistrySpec("DockerHashActorSpec") with FlatSpecLike with Matchers with BeforeAndAfterAll {
  behavior of "DockerRegistryActor"

  override protected lazy val registryFlows = List(
    new DockerHubRegistry(DockerRegistryConfig.default),
    new GcrRegistry(DockerRegistryConfig.default),
    new QuayRegistry(DockerRegistryConfig.default)
  )

  it should "retrieve a public docker hash" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("ubuntu:latest")
    
    expectMsgPF(5 second) {
      case DockerInfoSuccessResponse(DockerInformation(DockerHashResult(alg, hash), _), _) => 
        alg shouldBe "sha256"
        hash should not be empty
    }
  }

  it should "retrieve a public docker hash on gcr" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("gcr.io/google-containers/alpine-with-bash:1.0")

    expectMsgPF(5 second) {
      case DockerInfoSuccessResponse(DockerInformation(DockerHashResult(alg, hash), _), _) =>
        alg shouldBe "sha256"
        hash should not be empty
    }
  }
  
  it should "send image not found message back if the image does not exist" taggedAs IntegrationTest in {
    val notFound = makeRequest("ubuntu:nonexistingtag")
    dockerActor ! notFound

    expectMsgClass(5 seconds, classOf[DockerInfoNotFound])
  }

  it should "send unauthorized message back if the user doesn't have permission on this image" taggedAs IntegrationTest in {
    val unauthorized = makeRequest("tjeandet/sinatra:v1")
    dockerActor ! unauthorized

    expectMsgClass(5 seconds, classOf[DockerInfoUnauthorized])
  }
  
  it should "send an unrecognized host message if no flow can process the docker string" taggedAs IntegrationTest in {
    val unauthorized = makeRequest("unknown.io/image:v1")
    dockerActor ! unauthorized

    expectMsgClass(5 seconds, classOf[DockerHashUnknownRegistry])
  }
  
  it should "cache results" in {
    
    val image1 = dockerImage("ubuntu:latest")
    val request = DockerInfoRequest(image1)
    
    val hashSuccess = DockerHashResult("sha256", "hashvalue")
    val responseSuccess = DockerInfoSuccessResponse(DockerInformation(hashSuccess, None), request)
    val mockResponseSuccess = MockHashResponse(responseSuccess, 1)
    
    val responseFailure = DockerInfoFailedResponse(new Exception("Docker hash failed - part of test flow"), request)
    val mockResponseFailure = MockHashResponse(responseFailure, 1)
    
    // Send back success, failure, success, failure, ...
    val mockHttpFlow = new DockerRegistryMock(mockResponseSuccess, mockResponseFailure)
    val dockerActorWithCache = system.actorOf(DockerInfoActor.props(Seq(mockHttpFlow), 1000, 3 seconds, 10))
    
    dockerActorWithCache ! request
    expectMsg(DockerInfoSuccessResponse(DockerInformation(hashSuccess, None), request))
    // Necessary to give some time to the cache to be updated - as it's decoupled from sending back the response
    Thread.sleep(1000)
    
    dockerActorWithCache ! request
    // Without caching, the second request would have yielded a Failure since the mock flow alternates between a success and a failure
    // Getting a success here means the request didn't make it to the stream
    expectMsg(responseSuccess)
    // Second check to make sure only one of the requests made it to the stream
    mockHttpFlow.count() shouldBe 1
    // Wait 3 seconds to give the cache time to expire the value
    Thread.sleep(3000)
    dockerActorWithCache ! request
    // Expect a failure this time, as the cache should be empty
    expectMsg(responseFailure)
    mockHttpFlow.count() shouldBe 2
  }


  it should "not deadlock" taggedAs IntegrationTest in {
    lazy val dockerActorScale = system.actorOf(DockerInfoActor.props(registryFlows, 1000, 20.minutes, 0))
    0 until 400 foreach { _ =>
      dockerActorScale ! makeRequest("gcr.io/google-containers/alpine-with-bash:1.0")
    }

    val received = receiveN(400, 1 minute)
    received foreach { _ shouldBe a[DockerInfoSuccessResponse] }

    system.stop(dockerActorScale)
  }

  it should "get the manifest from a manifest list" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("ubuntu@sha256:6d0e0c26489e33f5a6f0020edface2727db9489744ecc9b4f50c7fa671f23c49")
    val response = expectMsgClass(classOf[DockerInfoSuccessResponse])
    response.dockerInformation.dockerCompressedSize shouldNot be(None)
  }
}
