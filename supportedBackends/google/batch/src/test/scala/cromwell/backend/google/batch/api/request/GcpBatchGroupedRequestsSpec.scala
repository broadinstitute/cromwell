package cromwell.backend.google.batch.api.request

import com.google.cloud.batch.v1.{CreateJobRequest, DeleteJobRequest, GetJobRequest}
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpec

class GcpBatchGroupedRequestsSpec extends AnyWordSpec with Matchers {

  "enqueue" should {
    "return a new instance" in {
      val group = GcpBatchGroupedRequests.empty
      val (newGroup, _) = group.enqueue(runRequest)
      newGroup shouldNot be(group)
    }

    "accept the possible requests" in {
      val _ = createGroup(runRequest, abortRequest, getRequest)

      succeed
    }
  }

  "size" should {
    "return the correct size" in {
      val group = createGroup(runRequest, runRequest, abortRequest, abortRequest, getRequest, getRequest)
      group.size shouldBe 6
    }
  }

  "entries" should {
    "return the entries in the insertion order" in {
      val a = runRequest
      val b = runRequest
      val group = createGroup(a, b)
      group.entries.map(_._1) shouldBe List(a, b)
    }
  }

  private def createGroup(requests: BatchApiRequest*): GcpBatchGroupedRequests =
    requests.foldLeft(GcpBatchGroupedRequests.empty) { case (acc, cur) =>
      acc.enqueue(cur)._1
    }

  private def runRequest =
    BatchRunCreationRequest(
      workflowId = WorkflowId.randomId(),
      requester = null,
      httpRequest = CreateJobRequest.newBuilder().build()
    )

  private def getRequest =
    BatchStatusPollRequest(
      workflowId = WorkflowId.randomId(),
      requester = null,
      httpRequest = GetJobRequest.newBuilder().build(),
      jobId = StandardAsyncJob(java.util.UUID.randomUUID().toString)
    )

  private def abortRequest =
    BatchAbortRequest(
      workflowId = WorkflowId.randomId(),
      requester = null,
      httpRequest = DeleteJobRequest.newBuilder().build(),
      jobId = StandardAsyncJob(java.util.UUID.randomUUID().toString)
    )
}
