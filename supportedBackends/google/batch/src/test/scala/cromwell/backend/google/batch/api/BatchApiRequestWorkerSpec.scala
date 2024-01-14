package cromwell.backend.google.batch.api

import akka.actor.{ActorRef, Props}
import akka.testkit._
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.{HttpHeaders, HttpRequest}
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchApiException,
  BatchApiRequestFailed,
  BatchStatusPollRequest,
  BatchWorkerRequestWork
}
import cromwell.backend.google.batch.api.TestBatchApiRequestWorker.{
  BatchApiBatchCallbackResponse,
  CallbackFailure,
  CallbackSuccess
}
import cromwell.backend.google.batch.models.{GcpBatchConfiguration, RunStatus}
import cromwell.core.{ExecutionEvent, TestKitSuite}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
//import io.grpc.Status
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.collection.immutable.Queue
import scala.concurrent.duration._
import scala.concurrent.{Future, Promise}
import scala.util.Try

class BatchApiRequestWorkerSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with Eventually
    with BeforeAndAfter {

  // TODO: Shall this be a string?
  implicit val batchHandler: TestBatchApiBatchHandler[String] = new TestBatchApiBatchHandler[String] {
    override def statusPollResultHandler(pollRequest: BatchStatusPollRequest,
                                         completionPromise: Promise[Try[Unit]]
    ): JsonBatchCallback[String] = new JsonBatchCallback[String] {
      override def onFailure(e: GoogleJsonError, responseHeaders: HttpHeaders): Unit =
        println(s"JsonBatchCallback -> onFailure: $e, responseHeaders=$responseHeaders")

      override def onSuccess(t: String, responseHeaders: HttpHeaders): Unit =
        println(s"JsonBatchCallback -> onSuccess: $t, responseHeaders=$responseHeaders")
    }
  }

  behavior of "BatchApiRequestWorker"

  implicit val TestExecutionTimeout: FiniteDuration = 10.seconds.dilated
  implicit val DefaultPatienceConfig: PatienceConfig = PatienceConfig(TestExecutionTimeout)
  val AwaitAlmostNothing: FiniteDuration = 30.milliseconds.dilated

  val managerProbe: TestProbe = TestProbe()
  var workerActor: TestActorRef[TestBatchApiRequestWorker] = _
//  TestActorRef(
//    TestBatchApiRequestWorker.props(managerProbe, config, registry)
//  )

  it should "correctly calculate batch intervals" in {
    import eu.timepit.refined.auto._
    BatchApiRequestManager.determineBatchInterval(10) should be(11111.milliseconds)
    BatchApiRequestManager.determineBatchInterval(100000) shouldBe 1.millisecond
  }

  it should "query for work and wait for a reply" in {
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[BatchApiRequestManager.BatchWorkerRequestWork])
    managerProbe.expectNoMessage(max = AwaitAlmostNothing)
  }

  it should "respond correctly with various run statuses" in {
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[BatchApiRequestManager.BatchWorkerRequestWork])

    val requester1 = TestProbe("requester1")
    val query1 = BatchStatusPollRequest(null, requester1.ref, null, null)
    val requester2 = TestProbe("requester2")
    val query2 = BatchStatusPollRequest(null, requester2.ref, null, null)
    val requester3 = TestProbe("requester3")
    val query3 = BatchStatusPollRequest(null, requester3.ref, null, null)

    // For two requests the callback succeeds (first with RunStatus.Success, then RunStatus.Failed). The third callback fails (simulating a network timeout, for example):
    batchHandler.callbackResponses :+= CallbackSuccess
    batchHandler.callbackResponses :+= CallbackSuccess
    batchHandler.callbackResponses :+= CallbackFailure

    val successStatus = RunStatus.Succeeded(Seq.empty[ExecutionEvent])
//    val successStatus = RunStatus.Success(Seq.empty[ExecutionEvent], None, None, None)
    val failureStatus = RunStatus.Failed(
      eventList = Seq.empty[ExecutionEvent]
    )
//    val failureStatus = RunStatus.UnsuccessfulRunStatus(
//      errorCode = Status.UNIMPLEMENTED,
//      errorMessage = Option.empty[String],
//      eventList = Seq.empty[ExecutionEvent],
//      machineType = None,
//      zone = None,
//      instanceName = None,
//      wasPreemptible = false
//    )
    batchHandler.operationStatusResponses :+= successStatus
    batchHandler.operationStatusResponses :+= failureStatus

    workerActor.tell(msg = BatchApiRequestManager.BatchApiWorkBatch(NonEmptyList(query1, List(query2, query3))),
                     sender = managerProbe.ref
    )
    eventually(batchHandler.runBatchRequested should be(true))

    // The manager shouldn't have been asked for more work yet:
    managerProbe.expectNoMessage(max = AwaitAlmostNothing)

    // Ok, let's trigger the callbacks:
    batchHandler.executeBatch()

    requester1.expectMsg(successStatus)
    requester2.expectMsg(failureStatus)
    requester3.expectNoMessage(max = AwaitAlmostNothing)

    // Requester3 expected nothing... Instead, the manager expects an API failure notification and then a request for more work:
    managerProbe.expectMsgPF(TestExecutionTimeout) { case failure: BatchApiRequestFailed =>
      if (!failure.cause.isInstanceOf[BatchApiException])
        fail("Unexpected failure cause class: " + failure.cause.getClass.getSimpleName)
      if (failure.query != query2 && failure.query != query3) fail("Unexpected query caused failure: " + failure.query)
    }
    managerProbe.expectMsg(BatchWorkerRequestWork(BatchApiRequestWorker.MaxBatchSize))
    managerProbe.expectNoMessage(max = AwaitAlmostNothing)
  }
}

//noinspection ScalaUnusedSymbol
class TestBatchApiRequestWorker(manager: ActorRef, qps: Int Refined Positive, registryProbe: ActorRef)(implicit
  batchHandler: TestBatchApiBatchHandler[_]
) extends BatchApiRequestWorker(manager, 10.milliseconds, registryProbe) {
  override def createBatch(): BatchRequest = null
  override def runBatch(batch: BatchRequest): Unit = batchHandler.runBatch()
}

//noinspection ScalaUnusedSymbol
abstract class TestBatchApiBatchHandler[O >: Null] extends GcpBatchApiRequestHandler {
  var operationStatusResponses: Queue[RunStatus] = Queue.empty
  var resultHandlers: Queue[JsonBatchCallback[O]] = Queue.empty
  var callbackResponses: Queue[BatchApiBatchCallbackResponse] = Queue.empty
  var runBatchRequested: Boolean = false

  def createBatch(genomicsInterface: O): BatchRequest = null
  def runBatch(): Unit = runBatchRequested = true

  def executeBatch(): Unit =
    resultHandlers.zip(callbackResponses) foreach { case (handler, response) =>
      response match {
        case CallbackSuccess => handler.onSuccess(null, null)
        case CallbackFailure =>
          val error: GoogleJsonError = new GoogleJsonError()
          handler.onFailure(error, null)
      }
    }

  def enqueueStatusPollInBatch(pollingRequest: BatchStatusPollRequest, batch: BatchRequest): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = statusPollResultHandler(pollingRequest, completionPromise)
    addStatusPollToBatch(null, batch, resultHandler)
    completionPromise.future
  }

  def statusPollResultHandler(pollRequest: BatchStatusPollRequest,
                              completionPromise: Promise[Try[Unit]]
  ): JsonBatchCallback[O]

  def addStatusPollToBatch(httpRequest: HttpRequest, batch: BatchRequest, resultHandler: JsonBatchCallback[O]): Unit =
    resultHandlers :+= resultHandler

  def mockStatusInterpreter(operation: O): RunStatus = {
    val (status, newQueue) = operationStatusResponses.dequeue
    operationStatusResponses = newQueue
    status
  }
}

object TestBatchApiRequestWorker {
  def props(manager: ActorRef, jesConfiguration: GcpBatchConfiguration, registryProbe: ActorRef)(implicit
    batchHandler: TestBatchApiBatchHandler[_]
  ): Props =
    Props(
      new TestBatchApiRequestWorker(
        manager,
        jesConfiguration.batchAttributes.qps,
        registryProbe
      )
    )

  sealed trait BatchApiBatchCallbackResponse
  case object CallbackSuccess extends BatchApiBatchCallbackResponse
  case object CallbackFailure extends BatchApiBatchCallbackResponse
}
