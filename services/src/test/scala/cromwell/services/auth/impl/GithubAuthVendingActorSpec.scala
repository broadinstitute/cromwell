package cromwell.services.auth.impl

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.TestProbe
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.services.auth.GithubAuthVending.{
  GithubAuthRequest,
  GithubAuthTokenResponse,
  GithubAuthVendingFailure,
  NoGithubAuthResponse
}
import cromwell.services.auth.ecm.EcmService
import cromwell.services.auth.impl.GithubAuthVendingActorSpec.TestGithubAuthVendingActor
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.{ExecutionContext, Future}

class GithubAuthVendingActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Eventually {

  final private val githubAuthEnabledServiceConfig =
    ConfigFactory.parseString(s"""
                                 |enabled = true
                                 |auth.azure = true
                                 |ecm.base-url = "https://mock-ecm-url.org"
        """.stripMargin)

  behavior of "GithubAuthVendingActor"

  it should "return NoGithubAuthResponse if GithubAuthVending is disabled" in {
    val serviceConfig = ConfigFactory.parseString(s"""
                                                     |enabled = false
                                                     |auth.azure = false
      """.stripMargin)

    val serviceRegistryActor = TestProbe()
    val actor = system.actorOf(
      Props(new GithubAuthVendingActor(serviceConfig, ConfigFactory.parseString(""), serviceRegistryActor.ref))
    )

    eventually {
      serviceRegistryActor.send(actor, GithubAuthRequest("valid_user_token"))
      serviceRegistryActor.expectMsg(NoGithubAuthResponse)
    }
  }

  it should "return invalid configuration error if no ECM base url is found" in {
    val serviceConfig = ConfigFactory.parseString(s"""
                                                     |enabled = true
                                                     |auth.azure = true
      """.stripMargin)

    val serviceRegistryActor = TestProbe()
    val actor = system.actorOf(
      Props(new GithubAuthVendingActor(serviceConfig, ConfigFactory.parseString(""), serviceRegistryActor.ref))
    )

    eventually {
      serviceRegistryActor.send(actor, GithubAuthRequest("valid_user_token"))
      serviceRegistryActor.expectMsg(
        GithubAuthVendingFailure("Invalid configuration for service 'GithubAuthVending': missing 'ecm.base-url' value.")
      )
    }
  }

  it should "return Github token if found" in {
    val serviceRegistryActor = TestProbe()
    val actor = system.actorOf(
      Props(
        new TestGithubAuthVendingActor(githubAuthEnabledServiceConfig,
                                       ConfigFactory.parseString(""),
                                       serviceRegistryActor.ref
        )
      )
    )

    eventually {
      serviceRegistryActor.send(actor, GithubAuthRequest("valid_user_token"))
      serviceRegistryActor.expectMsg(GithubAuthTokenResponse("gha_token"))
    }
  }

  it should "return failure message if ECM service returns non-successful response" in {
    val serviceRegistryActor = TestProbe()
    val actor = system.actorOf(
      Props(
        new TestGithubAuthVendingActor(githubAuthEnabledServiceConfig,
                                       ConfigFactory.parseString(""),
                                       serviceRegistryActor.ref
        )
      )
    )

    eventually {
      serviceRegistryActor.send(actor, GithubAuthRequest("invalid_user_token"))
      serviceRegistryActor.expectMsg(GithubAuthVendingFailure("Exception thrown for testing purposes"))
    }
  }
}

class TestEcmService(baseUrl: String) extends EcmService(baseUrl) {

  override def getGithubAccessToken(
    userToken: String
  )(implicit actorSystem: ActorSystem, ec: ExecutionContext): Future[String] =
    userToken match {
      case "valid_user_token" => Future.successful("gha_token")
      case _ => Future.failed(new RuntimeException("Exception thrown for testing purposes"))
    }
}

object GithubAuthVendingActorSpec {

  class TestGithubAuthVendingActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
      extends GithubAuthVendingActor(serviceConfig, globalConfig, serviceRegistryActor) {
    override lazy val ecmServiceOpt: Option[EcmService] = Some(new TestEcmService("https://mock-ecm-url.org"))
  }
}
