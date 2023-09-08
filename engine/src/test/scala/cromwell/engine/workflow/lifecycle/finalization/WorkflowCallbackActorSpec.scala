package cromwell.engine.workflow.lifecycle.finalization

import akka.testkit._
import akka.http.scaladsl.client.RequestBuilding.Post
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.model.{HttpResponse, StatusCodes}
import akka.testkit.TestProbe
import common.mock.MockSugar
import cromwell.core.retry.SimpleExponentialBackoff
import org.mockito.Mockito._
import cromwell.core.{TestKitSuite, WorkflowFailed, WorkflowId, WorkflowSucceeded}
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackActor.PerformCallbackCommand
import cromwell.engine.workflow.lifecycle.finalization.WorkflowCallbackJsonSupport._
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.util.{GracefulShutdownHelper, WomMocks}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import wom.values.WomString

import java.net.URI
import java.time.Instant
import scala.concurrent.duration._
import scala.concurrent.Future

class WorkflowCallbackActorSpec
  extends TestKitSuite with AnyFlatSpecLike with Matchers with MockSugar {

  behavior of "WorkflowCallbackActor"

  private implicit val ec = system.dispatcher

  private val msgWait = 10.second.dilated
  private val awaitAlmostNothing = 1.second
  private val serviceRegistryActor = TestProbe("testServiceRegistryActor")
  private val deathWatch = TestProbe("deathWatch")
  private val mockUri = new URI("http://example.com")
  private val basicOutputs = WomMocks.mockOutputExpectations(List("foo" -> WomString("bar")).toMap)

  private val httpSuccess = Future.successful(HttpResponse.apply(StatusCodes.OK))
  private val httpFailure = Future.successful(HttpResponse.apply(StatusCodes.GatewayTimeout))

  private def metadataEvents(workflowId: WorkflowId, successful: Boolean) = (
    MetadataEvent(
      MetadataKey(workflowId, None, "workflowCallback", "successful"),
      MetadataValue(successful)
    ),
    MetadataEvent(
      MetadataKey(workflowId, None, "workflowCallback", "url"),
      MetadataValue(mockUri.toString)
    ),
    MetadataEvent(
      MetadataKey(workflowId, None, "workflowCallback", "timestamp"),
      MetadataValue(Instant.now())
    )
  )

  it should "send a command to the default URI and record the correct metadata" in {
    // Setup
    val workflowId = WorkflowId.randomId()
    val mockHttpClient = mock[CallbackHttpHandler]

    val expectedPostBody = CallbackMessage(
      workflowId.toString,
      WorkflowSucceeded.toString,
      Option(basicOutputs.outputs.map(entry => (entry._1.name, entry._2)))
    )
    val expectedRequest = Post(mockUri.toString, expectedPostBody)
    val (expectedResultMetadata, expectedUriMetadata, expectedTimestampMetadata) = metadataEvents(workflowId, true)

    when(mockHttpClient.sendRequest(expectedRequest)).thenReturn(httpSuccess)

    val props = WorkflowCallbackActor.props(
      serviceRegistryActor.ref,
      Option(mockUri),
      httpClient = mockHttpClient
    )
    val workflowCallbackActor = system.actorOf(props, "testWorkflowCallbackActorSuccess")

    // Do the thing
    val cmd = PerformCallbackCommand(
      workflowId = workflowId, uri = None, terminalState = WorkflowSucceeded, workflowOutputs = basicOutputs
    )
    workflowCallbackActor ! cmd

    // Check the result
    serviceRegistryActor.expectMsgPF(msgWait) {
      case PutMetadataAction(List(resultEvent, uriEvent, timestampEvent), _) =>
        resultEvent.key shouldBe expectedResultMetadata.key
        resultEvent.value shouldBe expectedResultMetadata.value
        uriEvent.key shouldBe expectedUriMetadata.key
        uriEvent.value shouldBe expectedUriMetadata.value
        timestampEvent.key shouldBe expectedTimestampMetadata.key
        // Not checking timestamp value because it won't match
      case _ =>
    }

    verify(mockHttpClient, times(1)).sendRequest(expectedRequest)

    // Shut it down
    deathWatch.watch(workflowCallbackActor)
    workflowCallbackActor ! GracefulShutdownHelper.ShutdownCommand
    deathWatch.expectTerminated(workflowCallbackActor, msgWait)

  }

  it should "correctly handle a callback that fails, then succeeds" in {
    // Setup
    val workflowId = WorkflowId.randomId()
    val mockHttpClient = mock[CallbackHttpHandler]

    val expectedPostBody = CallbackMessage(
      workflowId.toString,
      WorkflowSucceeded.toString,
      Option(basicOutputs.outputs.map(entry => (entry._1.name, entry._2)))
    )
    val expectedRequest = Post(mockUri.toString, expectedPostBody)

    val (expectedResultMetadata, expectedUriMetadata, expectedTimestampMetadata) = metadataEvents(workflowId, true)

    when(mockHttpClient.sendRequest(expectedRequest)).thenReturn(httpFailure, httpFailure, httpSuccess)

    val props = WorkflowCallbackActor.props(
      serviceRegistryActor.ref,
      Option(mockUri),
      httpClient = mockHttpClient
    )
    val workflowCallbackActor = system.actorOf(props, "testWorkflowCallbackActorFailThenSuccees")

    // Do the thing
    val cmd = PerformCallbackCommand(
      workflowId = workflowId, uri = None, terminalState = WorkflowSucceeded, workflowOutputs = basicOutputs
    )
    workflowCallbackActor ! cmd

    // Check the result
    serviceRegistryActor.expectMsgPF(msgWait) {
      case PutMetadataAction(List(resultEvent, uriEvent, timestampEvent), _) =>
        resultEvent.key shouldBe expectedResultMetadata.key
        resultEvent.value shouldBe expectedResultMetadata.value
        uriEvent.key shouldBe expectedUriMetadata.key
        uriEvent.value shouldBe expectedUriMetadata.value
        timestampEvent.key shouldBe expectedTimestampMetadata.key
      // Not checking timestamp value because it won't match
      case _ =>
    }

    verify(mockHttpClient, times(3)).sendRequest(expectedRequest)

    // Shut it down
    deathWatch.watch(workflowCallbackActor)
    workflowCallbackActor ! GracefulShutdownHelper.ShutdownCommand
    deathWatch.expectTerminated(workflowCallbackActor, msgWait)

  }

  it should "retry sending the callback for a failing workflow before failing" in {
    // Setup
    val workflowId = WorkflowId.randomId()
    val mockHttpClient = mock[CallbackHttpHandler]

    val expectedPostBody = CallbackMessage(
      workflowId.toString,
      WorkflowFailed.toString,
      None
    )
    val expectedRequest = Post(mockUri.toString, expectedPostBody)

    val (expectedResultMetadata, expectedUriMetadata, expectedTimestampMetadata) = metadataEvents(workflowId, false)

    when(mockHttpClient.sendRequest(expectedRequest)).thenReturn(httpFailure)

    val props = WorkflowCallbackActor.props(
      serviceRegistryActor.ref,
      None,
      retryBackoff = Option(SimpleExponentialBackoff(500.millis, 1.minute, 1.1)),
      maxRetries = Option(5),
      httpClient = mockHttpClient
    )
    val workflowCallbackActor = system.actorOf(props, "testWorkflowCallbackActorFail")

    // Do the thing
    val cmd = PerformCallbackCommand(
      workflowId = workflowId, uri = Option(mockUri.toString), terminalState = WorkflowFailed, workflowOutputs = basicOutputs
    )
    workflowCallbackActor ! cmd

    // Check the result
    serviceRegistryActor.expectMsgPF(msgWait) {
      case PutMetadataAction(List(resultEvent, uriEvent, timestampEvent), _) =>
        resultEvent.key shouldBe expectedResultMetadata.key
        resultEvent.value shouldBe expectedResultMetadata.value
        uriEvent.key shouldBe expectedUriMetadata.key
        uriEvent.value shouldBe expectedUriMetadata.value
        timestampEvent.key shouldBe expectedTimestampMetadata.key
      // Not checking timestamp value because it won't match
      case _ =>
    }

    verify(mockHttpClient, times(5)).sendRequest(expectedRequest)

    // Shut it down
    deathWatch.watch(workflowCallbackActor)
    workflowCallbackActor ! GracefulShutdownHelper.ShutdownCommand
    deathWatch.expectTerminated(workflowCallbackActor, msgWait)

  }

  it should "not send a callback if no URI is provided" in {
    // Setup
    val workflowId = WorkflowId.randomId()
    val mockHttpClient = mock[CallbackHttpHandler]

    val props = WorkflowCallbackActor.props(
      serviceRegistryActor.ref,
      None,
      httpClient = mockHttpClient
    )
    val workflowCallbackActor = system.actorOf(props, "testWorkflowCallbackActorNoUri")

    // Do the thing
    val cmd = PerformCallbackCommand(
      workflowId = workflowId, uri = None, terminalState = WorkflowSucceeded, workflowOutputs = basicOutputs
    )
    workflowCallbackActor ! cmd

    // Check the result
    serviceRegistryActor.expectNoMessage(awaitAlmostNothing)
    verifyNoInteractions(mockHttpClient)

    // Shut it down
    deathWatch.watch(workflowCallbackActor)
    workflowCallbackActor ! GracefulShutdownHelper.ShutdownCommand
    deathWatch.expectTerminated(workflowCallbackActor, msgWait)
  }

  it should "no nothing when given a bogus URI" in {
    // Setup
    val workflowId = WorkflowId.randomId()
    val mockHttpClient = mock[CallbackHttpHandler]

    val props = WorkflowCallbackActor.props(
      serviceRegistryActor.ref,
      None,
      httpClient = mockHttpClient
    )
    val workflowCallbackActor = system.actorOf(props, "testWorkflowCallbackActorBogusUri")

    // Do the thing
    val cmd = PerformCallbackCommand(
      workflowId = workflowId,
      uri = Option("this is not a very good URI, is it"),
      terminalState = WorkflowSucceeded,
      workflowOutputs = basicOutputs
    )
    workflowCallbackActor ! cmd

    // Check the result
    serviceRegistryActor.expectNoMessage(awaitAlmostNothing)
    verifyNoInteractions(mockHttpClient)

    // Shut it down
    deathWatch.watch(workflowCallbackActor)
    workflowCallbackActor ! GracefulShutdownHelper.ShutdownCommand
    deathWatch.expectTerminated(workflowCallbackActor, msgWait)
  }
}
