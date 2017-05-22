package cromwell.core.callcaching.docker

import akka.http.scaladsl.Http
import akka.stream.ActorMaterializer
import akka.testkit.ImplicitSender
import cromwell.core.Tags.IntegrationTest
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.docker.DockerHashActor._
import cromwell.core.callcaching.docker.registryv2.flows.HttpFlowWithRetry.ContextWithRequest
import cromwell.core.callcaching.docker.registryv2.flows.dockerhub.DockerHubFlow
import cromwell.core.callcaching.docker.registryv2.flows.gcr.GoogleFlow
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class DockerHashActorSpec extends TestKitSuite with FlatSpecLike with Matchers with ImplicitSender {

  behavior of "DockerRegistryActor"
  
  implicit val materializer = ActorMaterializer()
  implicit val ex = system.dispatcher
  implicit val scheduler = system.scheduler
  
  val httpPool = Http().superPool[ContextWithRequest[DockerHashContext]]()
  val registryFlows = Seq(new DockerHubFlow(httpPool), new GoogleFlow(httpPool, 50000))
  // Disable cache by setting a cache size of 0 - A separate test tests the cache
  val dockerActor = system.actorOf(DockerHashActor.props(registryFlows, 1000, 20 minutes, 0)(materializer))
  
  private def dockerImage(string: String) = DockerImageIdentifier.fromString(string).get.asInstanceOf[DockerImageIdentifierWithoutHash]
  
  private def makeRequest(string: String) = {
    DockerHashRequest(dockerImage(string))
  }

  override protected def afterAll() = {
    system.stop(dockerActor)
    Http().shutdownAllConnectionPools()
    materializer.shutdown()
    super.afterAll()
  }
  
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


  it should "not deadlock" taggedAs IntegrationTest in {
    lazy val dockerActorScale = system.actorOf(DockerHashActor.props(registryFlows, 1000, 20.minutes, 0)(materializer))
    0 until 400 foreach { _ =>
      dockerActorScale ! makeRequest("gcr.io/google-containers/alpine-with-bash:1.0")
    }

    val received = receiveN(400, 1 minute)
    received foreach { _ shouldBe a[DockerHashResponseSuccess] }
  }
  
}
