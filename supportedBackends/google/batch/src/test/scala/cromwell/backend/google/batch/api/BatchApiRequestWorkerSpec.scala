package cromwell.backend.google.batch.api

import akka.actor.{ActorRef, PoisonPill}
import akka.testkit._
import cats.data.NonEmptyList
import com.google.cloud.batch.v1.{CancelJobRequest, CreateJobRequest, GetJobRequest}
import cromwell.backend.google.batch.api.BatchApiRequestManager.{
  BatchAbortRequest,
  BatchRunCreationRequest,
  BatchStatusPollRequest
}
import cromwell.backend.google.batch.api.request.{BatchApiRequestHandler, BatchRequestExecutor, GcpBatchGroupedRequests}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.{TestKitSuite, WorkflowId}
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import java.util.concurrent.atomic.AtomicReference
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

class BatchApiRequestWorkerSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with Eventually
    with BeforeAndAfter {

  behavior of "BatchApiRequestWorker"

  implicit val TestExecutionTimeout: FiniteDuration = 10.seconds.dilated
  implicit val DefaultPatienceConfig: PatienceConfig = PatienceConfig(TestExecutionTimeout)

  it should "work as expected" in {
    val managerProbe: TestProbe = TestProbe()
    val registry: TestProbe = TestProbe()
    val batchInterval = 1.second

    // dummy implementation that only enqueues the requests, ignoring the outcome
    implicit val batchHandler: BatchApiRequestHandler = new BatchApiRequestHandler {
      override def enqueue[T <: BatchApiRequestManager.BatchApiRequest](
        request: T,
        batchRequest: GcpBatchGroupedRequests,
        pollingManager: ActorRef
      )(implicit ec: ExecutionContext): GcpBatchGroupedRequests = batchRequest.enqueue(request)._1
    }

    // dummy implementations that count the number of processed requests
    import scala.language.reflectiveCalls
    val requestExecutor = new BatchRequestExecutor {
      val requestCount = new AtomicReference[Int]
      override def execute(groupedRequests: GcpBatchGroupedRequests)(implicit
        ec: ExecutionContext
      ): Future[List[Try[Unit]]] = Future {
        requestCount.updateAndGet(x => x + groupedRequests.entries.size)
        groupedRequests.entries.map(_ => Try(()))
      }
    }

    val workerActor: TestActorRef[BatchApiRequestWorker] = TestActorRef(
      BatchApiRequestWorker.props(
        pollingManager = managerProbe.ref,
        batchInterval = 1.second,
        serviceRegistryActor = registry.ref,
        requestExecutor
      )
    )

    // the actor just started, a message must be sent after the batchInterval
    managerProbe.expectMsgClass(max = batchInterval * 2, c = classOf[BatchApiRequestManager.BatchWorkerRequestWork])

    // but, no more messages are sent
    managerProbe.expectNoMessage(max = batchInterval)

    // until NoWorkToDo is sent, which schedules a BatchWorkerRequestWork message again
    workerActor ! BatchApiRequestManager.NoWorkToDo
    managerProbe.expectMsgClass(max = batchInterval * 2, c = classOf[BatchApiRequestManager.BatchWorkerRequestWork])

    // receiving BatchApiWorkBatch executes the requests, then, schedules a BatchWorkerRequestWork message if the work succeeds
    val requests = NonEmptyList(runCreationRequest, List(statusPollRequest, batchAbortRequest, statusPollRequest))
    workerActor ! BatchApiRequestManager.BatchApiWorkBatch(requests)

    managerProbe.expectMsgClass(max = batchInterval * 2, c = classOf[BatchApiRequestManager.BatchWorkerRequestWork])
    assert(requestExecutor.requestCount.get == requests.size,
           s"There were ${requests.size} expected requests but ${requestExecutor.requestCount.get} found"
    )
    requestExecutor.requestCount.set(0)

    // the same worker can process more work
    workerActor ! BatchApiRequestManager.BatchApiWorkBatch(requests)
    managerProbe.expectMsgClass(max = batchInterval * 2, c = classOf[BatchApiRequestManager.BatchWorkerRequestWork])
    assert(requestExecutor.requestCount.get == requests.size,
           s"There were ${requests.size} expected requests but ${requestExecutor.requestCount.get} found"
    )

    workerActor ! PoisonPill
  }

  private def statusPollRequest =
    BatchStatusPollRequest(
      WorkflowId.randomId(),
      requester = null,
      GetJobRequest.newBuilder().build(),
      StandardAsyncJob("id")
    )

  private def runCreationRequest =
    BatchRunCreationRequest(WorkflowId.randomId(), requester = null, CreateJobRequest.newBuilder().build())

  private def batchAbortRequest =
    BatchAbortRequest(
      WorkflowId.randomId(),
      requester = null,
      CancelJobRequest.newBuilder().build(),
      StandardAsyncJob("id")
    )

}
