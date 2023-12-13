package org.broadinstitute.dsde.workbench.cromwell.consumer

import au.com.dius.pact.consumer.dsl.LambdaDsl.newJsonBody
import au.com.dius.pact.consumer.dsl._
import au.com.dius.pact.consumer.{ConsumerPactBuilder, PactTestExecutionContext}
import au.com.dius.pact.core.model.RequestResponsePact
import cats.effect.IO
import cromwell.engine.workflow.lifecycle.finalization.{CallbackMessage, WorkflowCallbackJsonSupport}
import org.broadinstitute.dsde.workbench.cromwell.consumer.PactHelper._
import org.http4s.blaze.client.BlazeClientBuilder
import org.http4s.client.Client
import org.http4s.headers.Authorization
import org.http4s.{AuthScheme, Credentials, Uri}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import pact4s.scalatest.RequestResponsePactForger
import wom.values._

import java.util.concurrent.Executors
import scala.concurrent.ExecutionContext

class CbasClientSpec extends AnyFlatSpec with Matchers with RequestResponsePactForger {
  val ec: ExecutionContext = ExecutionContext.fromExecutorService(Executors.newCachedThreadPool())
  implicit val cs = IO.contextShift(ec)
  /*
    Define the folder that the pact contracts get written to upon completion of this test suite.
   */
  override val pactTestExecutionContext: PactTestExecutionContext =
    new PactTestExecutionContext(
      "./target/pacts"
    )

  val bearerToken = "my-token"
  val workflowId = "12345678-1234-1234-1111-111111111111"
  val completedState = "Succeeded"
  val workflowOutputs = Map(("wf.foo", WomString("bar")), ("wf.hello.hello", WomString("Hello")))
  val failures = List.empty[String]

  val workflowCallback: CallbackMessage = CallbackMessage(
    workflowId,
    completedState,
    workflowOutputs,
    failures
  )
  /*
      Get a json representation of a workflowCallback object.
   */
  val workflowCallbackJson = WorkflowCallbackJsonSupport.callbackMessageFormat.write(workflowCallback)

  val updateCompletedRunDsl = newJsonBody { o =>
    o.stringType("workflowId", workflowCallback.workflowId)
    o.stringType("state", workflowCallback.state)
    o.stringType("outputs", workflowCallbackJson.asJsObject.fields("outputs").toString())
    o.array("failures",
            f =>
              workflowCallback.failures.foreach(f.stringType)
    )
    ()
  }.build

  val consumerPactBuilder: ConsumerPactBuilder = ConsumerPactBuilder
    .consumer("cromwell")

  val pactProvider: PactDslWithProvider = consumerPactBuilder
    .hasPactWith("cbas")

  var pactUpdateCompletedRunDslResponse: PactDslResponse = buildInteraction(
    pactProvider,
    state = "post completed workflow results",
    uponReceiving = "Request to post workflow results",
    method = "POST",
    path = "/api/batch/v1/runs/results",
    requestHeaders = Seq("Authorization" -> "Bearer %s".format(bearerToken), "Content-type" -> "application/json"),
    requestBody = updateCompletedRunDsl,
    status = 200
  )
  override val pact: RequestResponsePact = pactUpdateCompletedRunDslResponse.toPact

  val client: Client[IO] =
    BlazeClientBuilder[IO](ExecutionContext.global).resource.allocated.unsafeRunSync()._1

  it should "successfully post workflow results" in {
    new CbasClientImpl[IO](client, Uri.unsafeFromString(mockServer.getUrl))
      .postWorkflowResults(
        Authorization(Credentials.Token(AuthScheme.Bearer, bearerToken)),
        CallbackMessage(workflowId, completedState, workflowOutputs, failures)
      )
      .attempt
      .unsafeRunSync() shouldBe Right(true)
  }
}
