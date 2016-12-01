package cromwell.backend.impl.jes.statuspolling

import akka.actor.{ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe}
import cromwell.core.{ExecutionEvent, TestKitSuite}
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}
import org.scalatest.concurrent.Eventually

import scala.concurrent.duration._
import akka.testkit._
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.{JesConfiguration, Run, RunStatus}
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.JesStatusPollQuery
import cromwell.backend.impl.jes.statuspolling.JesPollingActor.JesPollFailed
import cromwell.backend.impl.jes.statuspolling.TestJesPollingActor.{CallbackFailure, CallbackSuccess, JesBatchCallbackResponse}
import org.specs2.mock.Mockito

import scala.collection.immutable.Queue

class JesPollingActorSpec extends TestKitSuite("JesPollingActor") with FlatSpecLike with Matchers with Eventually with BeforeAndAfter with Mockito {

  behavior of "JesPollingActor"

  implicit val TestExecutionTimeout = 10.seconds.dilated
  implicit val DefaultPatienceConfig = PatienceConfig(TestExecutionTimeout)
  val AwaitAlmostNothing = 30.milliseconds.dilated

  import cromwell.backend.impl.jes.JesTestConfig.JesBackendConfigurationDescriptor
  val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

  var managerProbe: TestProbe = _
  var jpActor: TestActorRef[TestJesPollingActor] = _

  it should "correctly calculate batch intervals" in {
    JesPollingActor.determineBatchInterval(10) shouldBe 9.seconds
    JesPollingActor.determineBatchInterval(100) shouldBe 1.second
  }

  it should "query for work and wait for a reply" in {
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[JesApiQueryManager.RequestJesPollingWork])
    managerProbe.expectNoMsg(max = AwaitAlmostNothing)
  }

  it should "respond directly to requesters with various run statuses" in {
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[JesApiQueryManager.RequestJesPollingWork])

    val requester1 = TestProbe()
    val query1 = JesStatusPollQuery(requester1.ref, mock[Run])
    val requester2 = TestProbe()
    val query2 = JesStatusPollQuery(requester2.ref, mock[Run])
    val requester3 = TestProbe()
    val query3 = JesStatusPollQuery(requester3.ref, mock[Run])

    // For two requests the callback succeeds (first with RunStatus.Success, then RunStatus.Failed). The third callback fails:
    jpActor.underlyingActor.callbackResponses :+= CallbackSuccess
    jpActor.underlyingActor.callbackResponses :+= CallbackSuccess
    jpActor.underlyingActor.callbackResponses :+= CallbackFailure

    val successStatus = RunStatus.Success(Seq.empty[ExecutionEvent], None, None, None)
    val failureStatus = RunStatus.Failed(-1, List.empty[String], Seq.empty[ExecutionEvent], None, None, None)
    jpActor.underlyingActor.operationStatusResponses :+= successStatus
    jpActor.underlyingActor.operationStatusResponses :+= failureStatus

    jpActor.tell(msg = JesApiQueryManager.JesPollingWorkBatch(NonEmptyList(query1, List(query2, query3))), sender = managerProbe.ref)
    eventually { jpActor.underlyingActor.resultHandlers.size should be(3) }
    eventually { jpActor.underlyingActor.runBatchRequested should be(true) }

    // The manager shouldn't have been asked for more work yet:
    managerProbe.expectNoMsg(max = AwaitAlmostNothing)

    // Ok, let's trigger the callbacks:
    jpActor.underlyingActor.executeBatch()

    requester1.expectMsg(successStatus)
    requester2.expectMsg(failureStatus)
    requester3.expectMsgClass(classOf[JesPollFailed])

    // And the poller is done! Now the manager should now have (only one) request for more work:
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[JesApiQueryManager.RequestJesPollingWork])
  }

  before {
    managerProbe = TestProbe()
    jpActor = TestActorRef(TestJesPollingActor.props(managerProbe.ref, jesConfiguration), managerProbe.ref)
  }
}

object JesPollingActorSpec extends Mockito {
  def mockRun(runId: String): Run = {
    val run = mock[Run]
    run.runId returns runId
    run
  }
}

/**
  * Testable JES polling actor.
  * - Mocks out the methods which actually call out to JES, and allows the callbacks to be triggered in a testable way
  * - Also waits a **lot** less time before polls!
  */
class TestJesPollingActor(manager: ActorRef, qps: Int) extends JesPollingActor(manager, qps) with Mockito {

  override lazy val batchInterval = 10.milliseconds

  var operationStatusResponses: Queue[RunStatus] = Queue.empty
  var resultHandlers: Queue[JsonBatchCallback[Operation]] = Queue.empty
  var callbackResponses: Queue[JesBatchCallbackResponse] = Queue.empty
  var runBatchRequested: Boolean = false

  override private[statuspolling] def createBatch(run: Run): BatchRequest = null
  override private[statuspolling] def runBatch(batch: BatchRequest): Unit = runBatchRequested = true

  def executeBatch(): Unit = {
    resultHandlers.zip(callbackResponses) foreach { case (handler, response) => response match {
      case CallbackSuccess => handler.onSuccess(null, null)
      case CallbackFailure =>
        val error: GoogleJsonError = null
        handler.onFailure(error, null)
    }}
  }
  override private[statuspolling] def enqueueStatusPollInBatch(run: Run, batch: BatchRequest, resultHandler: JsonBatchCallback[Operation]): Unit = resultHandlers :+= resultHandler
  override private[statuspolling] def interpretOperationStatus(operation: Operation): RunStatus = {
    val (status, newQueue) = operationStatusResponses.dequeue
    operationStatusResponses = newQueue
    status
  }
  override private[statuspolling] def mkErrorString(e: GoogleJsonError) = "NA"
}

object TestJesPollingActor {
  def props(manager: ActorRef, jesConfiguration: JesConfiguration) = Props(new TestJesPollingActor(manager, jesConfiguration.qps))

  sealed trait JesBatchCallbackResponse
  case object CallbackSuccess extends JesBatchCallbackResponse
  case object CallbackFailure extends JesBatchCallbackResponse
}
