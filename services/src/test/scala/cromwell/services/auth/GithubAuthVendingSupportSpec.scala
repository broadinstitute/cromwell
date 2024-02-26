package cromwell.services.auth

import akka.actor.ActorRef
import akka.testkit.TestProbe
import akka.util.Timeout
import cromwell.core.TestKitSuite
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services.auth.GithubAuthVending.GithubAuthRequest
import cromwell.services.auth.GithubAuthVendingSupportSpec.TestGithubAuthVendingSupport
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.{Await, ExecutionContext, Future}
import scala.concurrent.duration._

class GithubAuthVendingSupportSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Eventually {

  behavior of "GithubAuthVendingSupport"

  it should "send a message to the service registry and handle success response" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)
    val provider = testSupport.importAuthProvider("user-token")
    val authHeader: Future[Map[String, String]] = provider.authHeader()

    serviceRegistryActor.expectMsg(GithubAuthRequest("user-token"))
    serviceRegistryActor.reply(GithubAuthVending.GithubAuthTokenResponse("github-token"))

    Await.result(authHeader, 10.seconds) should be(Map("Authorization" -> "Bearer github-token"))
  }

  it should "handle 'no auth' responses" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)
    val provider = testSupport.importAuthProvider("user-token")
    val authHeader: Future[Map[String, String]] = provider.authHeader()

    serviceRegistryActor.expectMsg(GithubAuthRequest("user-token"))
    serviceRegistryActor.reply(GithubAuthVending.NoGithubAuthResponse)

    Await.result(authHeader, 10.seconds) should be(Map.empty)
  }

  it should "convert error responses in Future failures" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)
    val provider = testSupport.importAuthProvider("user-token")
    val authHeader: Future[Map[String, String]] = provider.authHeader()

    serviceRegistryActor.expectMsg(GithubAuthRequest("user-token"))
    serviceRegistryActor.reply(GithubAuthVending.GithubAuthVendingFailure(new Exception("BOOM")))

    eventually {
      authHeader.isCompleted should be(true)
      authHeader.value.get.failed.get.getMessage should be("Failed to resolve github auth token")
      authHeader.value.get.failed.get.getCause.getMessage should be("BOOM")
    }
  }

  it should "handle timeouts" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref, 1.millisecond)
    val provider = testSupport.importAuthProvider("user-token")
    val authHeader: Future[Map[String, String]] = provider.authHeader()

    eventually {
      authHeader.isCompleted should be(true)
      authHeader.value.get.failed.get.getMessage should be("Unable to resolve github auth token within allowed time")
      authHeader.value.get.failed.get.getCause.getMessage should startWith("Ask timed out on")
    }
  }

  it should "handle ask failures" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)
    val provider = testSupport.importAuthProvider("user-token")
    val authHeader: Future[Map[String, String]] = provider.authHeader()

    serviceRegistryActor.expectMsg(GithubAuthRequest("user-token"))
    serviceRegistryActor.reply(ServiceRegistryFailure("GithubAuthVending"))

    eventually {
      authHeader.isCompleted should be(true)
      authHeader.value.get.failed.get.getMessage should be("Failed to resolve github auth token")
      // not the prettiest error message, but at least it should give us something to work with at debug time:
      authHeader.value.get.failed.get.getCause.getMessage should startWith("Cannot cast")
    }
  }

}

object GithubAuthVendingSupportSpec {
  class TestGithubAuthVendingSupport(val serviceRegistryActor: ActorRef, val timeout: Timeout = 10.seconds)
      extends GithubAuthVendingSupport {
    implicit override val ec: ExecutionContext = ExecutionContext.global
  }

}
