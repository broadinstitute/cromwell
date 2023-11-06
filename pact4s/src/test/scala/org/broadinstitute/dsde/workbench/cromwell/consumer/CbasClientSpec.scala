package org.broadinstitute.dsde.workbench.cromwell.consumer

import au.com.dius.pact.consumer.dsl.LambdaDsl.newJsonBody
import au.com.dius.pact.consumer.dsl._
import au.com.dius.pact.consumer.{ConsumerPactBuilder, PactTestExecutionContext}
import au.com.dius.pact.core.model.RequestResponsePact
import cats.effect.IO
import org.broadinstitute.dsde.workbench.cromwell.consumer.PactHelper._
import org.http4s.Uri
import org.http4s.blaze.client.BlazeClientBuilder
import org.http4s.client.Client
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import pact4s.scalatest.RequestResponsePactForger

import java.util.UUID
import java.util.concurrent.Executors
import scala.concurrent.ExecutionContext

class CbasClientSpec extends AnyFlatSpec with Matchers with RequestResponsePactForger {
  val ec = ExecutionContext.fromExecutorService(Executors.newCachedThreadPool())
  implicit val cs = IO.contextShift(ec)
  /*
    Define the folder that the pact contracts get written to upon completion of this test suite.
   */
  override val pactTestExecutionContext: PactTestExecutionContext =
    new PactTestExecutionContext(
      "./target/pacts"
    )

  val bearerToken = "my-token"
  val workflowId = "33333333-3333-3333-3333-333333333333"
  val workspaceId = "44444444-3333-3333-3333-333333333333"

  val state = "Aborted"
  val outputs = "[]"
  val failures = List.empty[String]

  val updateRequestDsl = newJsonBody { o =>
    o.stringType("workflowId", workflowId)
    o.stringType("state", state)
    o.stringType("outputs", outputs)
    o.array("failures",
      { f =>
        failures.foreach(f.stringType)
      })
    ()
  }.build

  val consumerPactBuilder: ConsumerPactBuilder = ConsumerPactBuilder
    .consumer("cromwell")

  val pactProvider: PactDslWithProvider = consumerPactBuilder
    .hasPactWith("cbas")


  var pactDslResponse: PactDslResponse = buildInteraction(
    pactProvider,
    uponReceiving = "Request to post workflow results",
    method = "POST",
    path = "/api/batch/v1/runs/results",
    requestHeaders = Seq("Authorization" -> "Bearer %s".formatted(bearerToken), "Content-type" -> "application/json"),
    requestBody = updateRequestDsl,
    status = 200
  )
  override val pact: RequestResponsePact = pactDslResponse.toPact

  val client: Client[IO] = {
    BlazeClientBuilder[IO](ExecutionContext.global).resource.allocated.unsafeRunSync()._1
  }

  it should "post workflow results" in {
    new CbasClientImpl[IO](client, Uri.unsafeFromString(mockServer.getUrl))
      .postWorkflowResults(bearerToken, UUID.fromString(workflowId), state, outputs, failures)
      .attempt
      .unsafeRunSync() shouldBe Right(true)
  }
}
