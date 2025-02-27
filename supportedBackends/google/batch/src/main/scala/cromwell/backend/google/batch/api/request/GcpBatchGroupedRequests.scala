package cromwell.backend.google.batch.api.request

import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.google.batch.api.BatchApiResponse

import scala.concurrent.{Future, Promise}
import scala.util.Try

// Mirrors com.google.api.client.googleapis.batch.BatchRequest but this is immutable
class GcpBatchGroupedRequests(requests: List[(BatchApiRequest, Promise[Try[BatchApiResponse]])]) {

  def enqueue(that: BatchApiRequest): (GcpBatchGroupedRequests, Future[Try[BatchApiResponse]]) = {
    val promise = Promise[Try[BatchApiResponse]]()
    val newRequests = (that, promise) :: requests
    new GcpBatchGroupedRequests(newRequests) -> promise.future
  }

  def size: Int = requests.size
  def entries: List[(BatchApiRequest, Promise[Try[BatchApiResponse]])] = requests.reverse

}

object GcpBatchGroupedRequests {
  def empty: GcpBatchGroupedRequests = new GcpBatchGroupedRequests(List.empty)
}
