package org.broadinstitute.dsde.workbench.cromwell.consumer

import akka.testkit._
import au.com.dius.pact.consumer.dsl._
import au.com.dius.pact.consumer.{ConsumerPactBuilder, PactTestExecutionContext}
import au.com.dius.pact.core.model.RequestResponsePact
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{CallOutputs, TestKitSuite, WorkflowId, WorkflowSucceeded}
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackActor.PerformCallbackCommand
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackConfig.StaticTokenAuth
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackJsonSupport._
import cromwell.engine.workflow.lifecycle.finalization.{CallbackMessage, WorkflowCallbackActor, WorkflowCallbackConfig}
import cromwell.services.metadata.MetadataKey
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.util.GracefulShutdownHelper
import org.broadinstitute.dsde.workbench.cromwell.consumer.PactHelper._
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import pact4s.scalatest.RequestResponsePactForger
import pact4s.sprayjson.implicits._
import wom.graph.GraphNodePort.GraphNodeOutputPort
import wom.graph.WomIdentifier
import wom.types.WomStringType
import wom.values._

import java.net.URI
import java.util.UUID
import scala.concurrent.ExecutionContextExecutor
import scala.concurrent.duration._

class CbasCallbackSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with RequestResponsePactForger {

  implicit val ec: ExecutionContextExecutor = system.dispatcher

  // Akka test setup
  private val msgWait = 10.second.dilated
  private val serviceRegistryActor = TestProbe("testServiceRegistryActor")
  private val deathWatch = TestProbe("deathWatch")

  private val callbackEndpoint = "/api/batch/v1/runs/results"
  private val bearerToken = "my-token"
  private val workflowId = WorkflowId(UUID.fromString("12345678-1234-1234-1111-111111111111"))
  private val basicOutputs = CallOutputs(
    Map(
      GraphNodeOutputPort(WomIdentifier("foo", "wf.foo"), WomStringType, null) -> WomString("bar"),
      GraphNodeOutputPort(WomIdentifier("hello", "wf.hello.hello"), WomStringType, null) -> WomString("Hello")
    )
  )

  // This is the message that we expect Cromwell to send CBAS in this test
  private val expectedCallbackMessage = CallbackMessage(
    workflowId.toString,
    "Succeeded",
    Map(("wf.foo", WomString("bar")), ("wf.hello.hello", WomString("Hello"))),
    List.empty
  )

  // Define the folder that the pact contracts get written to upon completion of this test suite.
  override val pactTestExecutionContext: PactTestExecutionContext =
    new PactTestExecutionContext(
      "./target/pacts"
    )

  val consumerPactBuilder: ConsumerPactBuilder = ConsumerPactBuilder
    .consumer("cromwell")

  val pactProvider: PactDslWithProvider = consumerPactBuilder
    .hasPactWith("cbas")

  var pactUpdateCompletedRunDslResponse: PactDslResponse = buildInteraction(
    pactProvider,
    state = "post completed workflow results",
    uponReceiving = "Request to post workflow results",
    method = "POST",
    path = callbackEndpoint,
    requestHeaders = Seq("Authorization" -> "Bearer %s".format(bearerToken), "Content-type" -> "application/json"),
    requestBody = expectedCallbackMessage,
    status = 200
  )
  override val pact: RequestResponsePact = pactUpdateCompletedRunDslResponse.toPact

  it should "send the right callback to the right URI" in {
    // Create actor
    val callbackConfig = WorkflowCallbackConfig.empty
      .copy(enabled = true)
      .copy(retryBackoff = SimpleExponentialBackoff(100.millis, 200.millis, 1.1))
      .copy(authMethod = Option(StaticTokenAuth(bearerToken)))
      .copy(defaultUri = Option(new URI(mockServer.getUrl + callbackEndpoint)))

    val props = WorkflowCallbackActor.props(
      serviceRegistryActor.ref,
      callbackConfig
    )
    val workflowCallbackActor = system.actorOf(props, "testWorkflowCallbackActorPact")

    // Send a command to trigger callback
    val cmd = PerformCallbackCommand(
      workflowId = workflowId,
      uri = None,
      terminalState = WorkflowSucceeded,
      workflowOutputs = basicOutputs,
      List.empty
    )
    workflowCallbackActor ! cmd

    // Confirm the callback was successful
    serviceRegistryActor.expectMsgPF(msgWait) {
      case PutMetadataAction(List(resultEvent, _, _), _) =>
        resultEvent.key shouldBe MetadataKey(workflowId, None, "workflowCallback", "successful")
      case _ =>
    }

    // Shut the actor down
    deathWatch.watch(workflowCallbackActor)
    workflowCallbackActor ! GracefulShutdownHelper.ShutdownCommand
    deathWatch.expectTerminated(workflowCallbackActor, msgWait)
  }
}
