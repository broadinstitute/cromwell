package cromwell.services.auth

import akka.actor.ActorRef
import akka.testkit.TestProbe
import akka.util.Timeout
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import cromwell.languages.util.ImportResolver.GithubImportAuthProvider
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services.auth.GithubAuthVending.{GithubAuthRequest, GithubToken, TerraToken}
import cromwell.services.auth.GithubAuthVendingSupportSpec.TestGithubAuthVendingSupport
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.{Await, ExecutionContext, Future}
import scala.concurrent.duration._

class GithubAuthVendingSupportSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Eventually {

  private def azureGithubAuthVendingConfig(enabled: Boolean = true) = ConfigFactory
    .parseString(
      s"""
         |services {
         |  GithubAuthVending {
         |    config {
         |      auth.azure = ${enabled}
         |    }
         |  }
         |}
         |""".stripMargin
    )

  implicit val timeout = Timeout(10.seconds)

  behavior of "GithubAuthVendingSupport"

  it should "send a message to the service registry and handle success response" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)
    val provider = testSupport.importAuthProvider("user-token")
    val authHeader: Future[Map[String, String]] = provider.authHeader()

    serviceRegistryActor.expectMsg(GithubAuthRequest(TerraToken("user-token")))
    serviceRegistryActor.reply(GithubAuthVending.GithubAuthTokenResponse(GithubToken("github-token")))

    Await.result(authHeader, 10.seconds) should be(Map("Authorization" -> "Bearer github-token"))
  }

  it should "handle 'no auth' responses" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)
    val provider = testSupport.importAuthProvider("user-token")
    val authHeader: Future[Map[String, String]] = provider.authHeader()

    serviceRegistryActor.expectMsg(GithubAuthRequest(TerraToken("user-token")))
    serviceRegistryActor.reply(GithubAuthVending.NoGithubAuthResponse)

    Await.result(authHeader, 10.seconds) should be(Map.empty)
  }

  it should "convert error responses in Future failures" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)
    val provider = testSupport.importAuthProvider("user-token")
    val authHeader: Future[Map[String, String]] = provider.authHeader()

    serviceRegistryActor.expectMsg(GithubAuthRequest(TerraToken("user-token")))
    serviceRegistryActor.reply(GithubAuthVending.GithubAuthVendingFailure("BOOM"))

    eventually {
      authHeader.isCompleted should be(true)
      authHeader.value.get.failed.get.getMessage should be("Failed to resolve GitHub auth token. Error: BOOM")
    }
  }

  it should "handle timeouts" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)
    val provider = testSupport.importAuthProvider("user-token")(Timeout(1.millisecond))
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

    serviceRegistryActor.expectMsg(GithubAuthRequest(TerraToken("user-token")))
    serviceRegistryActor.reply(ServiceRegistryFailure("GithubAuthVending"))

    eventually {
      authHeader.isCompleted should be(true)
      authHeader.value.get.failed.get.getMessage should be("Failed to resolve github auth token")
      // not the prettiest error message, but at least it should give us something to work with at debug time:
      authHeader.value.get.failed.get.getCause.getMessage should startWith("Cannot cast")
    }
  }

  it should "return Github import auth provider when Azure auth is enabled" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)

    testSupport.importAuthProvider(azureGithubAuthVendingConfig()) match {
      case Valid(providerOpt) =>
        providerOpt.isEmpty shouldBe false
        providerOpt.get.isInstanceOf[GithubImportAuthProvider] shouldBe true
        providerOpt.get.validHosts shouldBe List("github.com", "githubusercontent.com", "raw.githubusercontent.com")
      case Invalid(e) => fail(s"Unexpected failure: $e")
    }
  }

  it should "return no import auth provider when Azure auth is disabled" in {
    val serviceRegistryActor = TestProbe()
    val testSupport = new TestGithubAuthVendingSupport(serviceRegistryActor.ref)

    testSupport.importAuthProvider(azureGithubAuthVendingConfig(false)) match {
      case Valid(providerOpt) => providerOpt.isEmpty shouldBe true
      case Invalid(e) => fail(s"Unexpected failure: $e")
    }
  }
}

object GithubAuthVendingSupportSpec {
  class TestGithubAuthVendingSupport(val serviceRegistryActor: ActorRef) extends GithubAuthVendingSupport {
    implicit override val ec: ExecutionContext = ExecutionContext.global
  }

}
