package org.broadinstitute.dsde.workbench.cromwell.consumer

import au.com.dius.pact.consumer.dsl._
import au.com.dius.pact.consumer.{ConsumerPactBuilder, PactTestExecutionContext}
import au.com.dius.pact.core.model.RequestResponsePact
import cats.effect.IO
import org.broadinstitute.dsde.workbench.cromwell.consumer.AuthHelper._
import org.broadinstitute.dsde.workbench.cromwell.consumer.PactHelper._
import org.broadinstitute.dsde.workbench.model.WorkbenchEmail
import org.http4s.Uri
import org.http4s.blaze.client.BlazeClientBuilder
import org.http4s.client.Client
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import pact4s.scalatest.RequestResponsePactForger

import java.util.concurrent.Executors
import scala.concurrent.ExecutionContext

class DrsHubClientSpec extends AnyFlatSpec with Matchers with RequestResponsePactForger {
  val ec = ExecutionContext.fromExecutorService(Executors.newCachedThreadPool())
  implicit val cs = IO.contextShift(ec)
  /*
    Define the folder that the pact contracts get written to upon completion of this test suite.
   */
  override val pactTestExecutionContext: PactTestExecutionContext =
    new PactTestExecutionContext(
      "./target/pacts"
    )

  // Uncomment this so that mock server will run on specific port (e.g. 9003) instead of dynamically generated port.
  // override val mockProviderConfig: MockProviderConfig = MockProviderConfig.httpConfig("localhost", 9003)

  // --- End of fixtures section

  // ---- Dsl for specifying pacts between consumer and provider
  // Lambda Dsl: required for generating matching rules.
  // Favored over old-style Pact Dsl using PactDslJsonBody.
  // This rule expects DrsHub to respond with
  // 1. ok status
  // 2. ok statuses matching the given subsystem states
  val consumerPactBuilder: ConsumerPactBuilder = ConsumerPactBuilder
    .consumer("cromwell-consumer")

  val pactProvider: PactDslWithProvider = consumerPactBuilder
    .hasPactWith("drshub-provider")

  // stateParams provides the desired subsystem states
  // for drshub provider to generate the expected response
  var pactDslResponse: PactDslResponse = buildInteraction(
    pactProvider,
    state = "Drshub is ok",
    stateParams = Map(),
    uponReceiving = "Request to get Drshub ok status",
    method = "GET",
    path = "/status",
    requestHeaders = Seq("Accept" -> "application/json"),
    status = 200,
    responseHeaders = Seq("Content-type" -> "application/json")
  )

  override val pact: RequestResponsePact = pactDslResponse.toPact

  val client: Client[IO] =
    BlazeClientBuilder[IO](ExecutionContext.global).resource.allocated.unsafeRunSync()._1

  /*
  we should use these tests to ensure that our client class correctly handles responses from the provider - i.e. decoding, error mapping, validation
   */
  it should "get DrsHub ok status" in {
    new DrsHubClientImpl[IO](client, Uri.unsafeFromString(mockServer.getUrl), mockAuthToken(WorkbenchEmail("")))
      .fetchSystemStatus()
      .attempt
      .unsafeRunSync() shouldBe Right(true)
  }
}
