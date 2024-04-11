package org.broadinstitute.dsde.workbench.cromwell.consumer

import au.com.dius.pact.consumer.dsl._
import au.com.dius.pact.consumer.{ConsumerPactBuilder, PactTestExecutionContext}
import au.com.dius.pact.core.model.RequestResponsePact
import cromwell.core.TestKitSuite
import cromwell.services.auth.GithubAuthVending.{GithubToken, TerraToken}
import cromwell.services.auth.ecm.EcmService
import org.scalatest.flatspec.AnyFlatSpecLike

import scala.concurrent.{Await, ExecutionContextExecutor}
import org.scalatest.matchers.should.Matchers
import pact4s.scalatest.RequestResponsePactForger

import scala.concurrent.duration.DurationInt

class EcmServiceContractSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with RequestResponsePactForger {

  implicit val ec: ExecutionContextExecutor = system.dispatcher

  /*
    Define the folder that the pact contracts get written to upon completion of this test suite.
   */
  override val pactTestExecutionContext: PactTestExecutionContext =
    new PactTestExecutionContext(
      "./target/pacts"
    )

  val consumerPactBuilder: ConsumerPactBuilder = ConsumerPactBuilder
    .consumer("cromwell")

  val pactProvider: PactDslWithProvider = consumerPactBuilder
    .hasPactWith("ecm")

  var testUser: String = "cromwell_test_user@test.com"
  var testBearerToken: String = "cromwellBearerToken"
  var testGithubToken: String = "githubToken"

  var pactDslResponse: PactDslResponse = pactProvider
    .`given`("a user is registered",
             Map[String, String](
               "userEmail" -> testUser,
               "bearerToken" -> testBearerToken
             )
    )
    .uponReceiving("a github token request")
    .method("GET")
    .path("/api/oauth/v1/github/access-token")
    .headers(Map[String, String]("Authorization" -> s"Bearer $testBearerToken"))
    .willRespondWith()
    .status(200)
    .bodyMatchingContentType("text/plain", testGithubToken)

  override val pact: RequestResponsePact = pactDslResponse.toPact

  it should "get a github token" in {
    Await.result(
      new EcmService(mockServer.getUrl)
        .getGithubAccessToken(TerraToken(testBearerToken)),
      10.seconds
    ) shouldBe GithubToken(testGithubToken)
  }
}
