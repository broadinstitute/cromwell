package cromwell.docker

import cromwell.core.Tags.IntegrationTest
import cromwell.docker.DockerHashActor._
import cromwell.docker.registryv2.flows.dockerhub.DockerHubFlow
import cromwell.docker.registryv2.flows.gcr.GoogleFlow
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class DockerHashActorSpec extends DockerFlowSpec("DockerHashActorSpec") with FlatSpecLike with Matchers {
  behavior of "DockerRegistryActor"

  override protected lazy val registryFlows = Seq(new DockerHubFlow(httpPool), new GoogleFlow(httpPool, 50000))

  it should "retrieve a public docker hash" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("ubuntu:latest")
    
    expectMsgPF(5 second) {
      case DockerHashResponseSuccess(DockerHashResult(alg, hash), _) => 
        alg shouldBe "sha256"
        hash should not be empty
    }
  }

  it should "retrieve a public docker hash on gcr" taggedAs IntegrationTest in {
    dockerActor ! makeRequest("gcr.io/google-containers/alpine-with-bash:1.0")

    expectMsgPF(5 second) {
      case DockerHashResponseSuccess(DockerHashResult(alg, hash), _) =>
        alg shouldBe "sha256"
        hash should not be empty
    }
  }
  
  it should "send image not found message back if the image does not exist" taggedAs IntegrationTest in {
    val notFound = makeRequest("ubuntu:nonexistingtag")
    dockerActor ! notFound

    expectMsgClass(5 seconds, classOf[DockerHashNotFound])
  }

  it should "send unauthorized message back if the user doesn't have permission on this image" taggedAs IntegrationTest in {
    val unauthorized = makeRequest("tjeandet/sinatra:v1")
    dockerActor ! unauthorized

    expectMsgClass(5 seconds, classOf[DockerHashUnauthorized])
  }
  
  it should "send an unrecognized host message if no flow can process the docker string" taggedAs IntegrationTest in {
    val unauthorized = makeRequest("unknown.io/image:v1")
    dockerActor ! unauthorized

    expectMsgClass(5 seconds, classOf[DockerHashUnknownRegistry])
  }
  
  it should "cache results" in {
    
    val image1 = dockerImage("ubuntu:latest")
    val request = DockerHashRequest(image1)
    
    val hashSuccess = DockerHashResult("sha256", "hashvalue")
    val responseSuccess = DockerHashResponseSuccess(hashSuccess, request)
    val mockResponseSuccess = MockHashResponse(responseSuccess, 1)
    
    val responseFailure = DockerHashFailedResponse(new Exception("Docker hash failed - part of test flow"), request)
    val mockResponseFailure = MockHashResponse(responseFailure, 1)
    
    // Send back success, failure, success, failure, ...
    val mockHttpFlow = new DockerFlowMock(mockResponseSuccess, mockResponseFailure)
    val dockerActorWithCache = system.actorOf(DockerHashActor.props(Seq(mockHttpFlow), 1000, 3 seconds, 10)(materializer))
    
    dockerActorWithCache ! request
    expectMsg(DockerHashResponseSuccess(hashSuccess, request))
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
  
}
