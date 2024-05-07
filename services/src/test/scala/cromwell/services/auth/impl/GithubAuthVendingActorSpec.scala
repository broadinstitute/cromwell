package cromwell.services.auth.impl

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.TestProbe
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.services.auth.GithubAuthVending.{
  GithubAuthRequest,
  GithubAuthTokenResponse,
  GithubAuthVendingFailure,
  GithubToken,
  NoGithubAuthResponse,
  TerraToken
}
import cromwell.services.auth.ecm.EcmService
import cromwell.services.auth.impl.GithubAuthVendingActorSpec.TestGithubAuthVendingActor
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.concurrent.{ExecutionContext, Future}

class GithubAuthVendingActorSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with Eventually
    with TableDrivenPropertyChecks {

  final private val validUserToken = "valid_user_token"
  final private val githubAuthEnabledServiceConfig =
    ConfigFactory.parseString(s"""
                                 |enabled = true
                                 |auth.azure = true
                                 |ecm.base-url = "https://mock-ecm-url.org"
        """.stripMargin)

  private def githubAuthConfigWithoutEcmUrl(flag: Boolean): Config =
    ConfigFactory.parseString(s"""
                                 |enabled = $flag
                                 |auth.azure = $flag
        """.stripMargin)

  private val testCases = Table(
    ("test name", "service config", "use TestEcmService class", "terra token", "response"),
    ("return NoGithubAuthResponse if GithubAuthVending is disabled",
     githubAuthConfigWithoutEcmUrl(false),
     false,
     validUserToken,
     NoGithubAuthResponse
    ),
    ("return invalid configuration error if no ECM base url is found",
     githubAuthConfigWithoutEcmUrl(true),
     false,
     validUserToken,
     GithubAuthVendingFailure("Invalid configuration for service 'GithubAuthVending': missing 'ecm.base-url' value.")
    ),
    ("return Github token if found",
     githubAuthEnabledServiceConfig,
     true,
     validUserToken,
     GithubAuthTokenResponse(GithubToken("gha_token"))
    ),
    ("return failure message if ECM service returns non-successful response",
     githubAuthEnabledServiceConfig,
     true,
     "invalid_user_token",
     GithubAuthVendingFailure("Exception thrown for testing purposes")
    )
  )

  behavior of "GithubAuthVendingActor"

  forAll(testCases) { (testName, serviceConfig, useTestEcmService, userTerraToken, expectedResponseMsg) =>
    it should testName in {
      val serviceRegistryActor = TestProbe()
      val actor = system.actorOf(
        Props(new TestGithubAuthVendingActor(serviceConfig, serviceRegistryActor.ref, useTestEcmService))
      )

      serviceRegistryActor.send(actor, GithubAuthRequest(TerraToken(userTerraToken)))
      eventually {
        serviceRegistryActor.expectMsg(expectedResponseMsg)
      }
    }
  }
}

class TestEcmService(baseUrl: String) extends EcmService(baseUrl) {

  override def getGithubAccessToken(
    userToken: TerraToken
  )(implicit actorSystem: ActorSystem, ec: ExecutionContext): Future[GithubToken] =
    userToken.value match {
      case "valid_user_token" => Future.successful(GithubToken("gha_token"))
      case _ => Future.failed(new RuntimeException("Exception thrown for testing purposes"))
    }
}

object GithubAuthVendingActorSpec {

  class TestGithubAuthVendingActor(serviceConfig: Config,
                                   serviceRegistryActor: ActorRef,
                                   useTestEcmServiceClass: Boolean
  ) extends GithubAuthVendingActor(serviceConfig, ConfigFactory.parseString(""), serviceRegistryActor) {

    override lazy val ecmServiceOpt: Option[EcmService] =
      if (useTestEcmServiceClass) Some(new TestEcmService("https://mock-ecm-url.org"))
      else ecmConfigOpt.baseUrl.map(url => new EcmService(url))
  }
}
