//package cromwell.backend.google.batch.api
//
//import akka.actor.ActorRef
//import cromwell.backend.google.batch.api.request.GcpBatchGroupedRequests
////import com.google.api.client.googleapis.batch.BatchRequest
////import com.google.api.client.googleapis.batch.json.JsonBatchCallback
////import com.google.api.client.googleapis.json.GoogleJsonError
////import com.google.api.client.http.HttpHeaders
//import com.google.api.gax.rpc.FixedHeaderProvider
//import com.google.cloud.batch.v1._
//import com.google.common.collect.ImmutableMap
//import com.google.longrunning.Operation
////import common.util.Backoff
//import cromwell.backend.google.batch.api.BatchApiRequestManager._
////import cromwell.backend.standard.StandardAsyncJob
////import cromwell.core.WorkflowId
////import cromwell.core.retry.SimpleExponentialBackoff
//
////import scala.concurrent.duration.DurationInt
//import scala.concurrent.{ExecutionContext, Future, Promise}
//import scala.util.Try
//
//// A batch of requests to invoke GCP Batch service
////
//// GCP Batch does not support batching many requests, we are simulating a way to do it.
//// NOTE: This is not thread-safe
//// TODO: Alex - what about the performance? PAPIv2 is backed by an ArrayList
////case class GcpBatchGroupedRequests(private val requests: List[BatchApiRequest]) {
////  // prepend operation is more efficient
////  def +(that: BatchApiRequest): GcpBatchGroupedRequests = copy(requests = that :: requests)
////
////  // preserve insertion order
////  def getRequests: List[BatchApiRequest] = requests.reverse
////}
//
//class GcpBatchApiRequestHandler {
//  def submit(request: CreateJobRequest): Job = withClient { client =>
//    client.createJobCallable
//      .call(request)
//  }
//
//  // TODO: Similar to PipelinesApiRequestHandler, this is supposed to not execute the request just yet but enqueue it
//  def enqueue(request: BatchApiRequest, batchRequest: GcpBatchGroupedRequests, pollingManager: ActorRef)(implicit
//    ec: ExecutionContext
//  ): Future[Try[Unit]] = request match {
//    case r: BatchRunCreationRequest =>
//      val completionPromise = Promise[Try[Unit]]()
//      println(s"ec: ${ec.hashCode()}")
////      val resultHandler = runCreationResultHandler(r, completionPromise, pollingManager)
////      // TODO: Alex - implement this
////      println()
////      addRunCreationToBatch(r, batchRequest, resultHandler)
////      batchRequest.add(r)
//      println(s"enqueue -> BatchRunCreationRequest: ${r}")
//      completionPromise.future
//
//    case r: BatchStatusPollRequest =>
//      println(s"enqueue -> BatchStatusPollRequest: ${r}")
//      val completionPromise = Promise[Try[Unit]]()
//      completionPromise.future
//
//    case r: BatchAbortRequest =>
//      println(s"enqueue -> BatchAbortRequest: ${r}")
//      val completionPromise = Promise[Try[Unit]]()
//      completionPromise.future
//  }
//
////  private def runCreationResultHandler(request: BatchRunCreationRequest,
////                                       completionPromise: Promise[Try[Unit]],
////                                       pollingManager: ActorRef
////  ) = {}
//
//  def makeBatchRequest: GcpBatchGroupedRequests = new GcpBatchGroupedRequests
//
//  def query(request: GetJobRequest): Job = withClient { client =>
//    client.getJob(request)
//  }
//
//  def abort(request: DeleteJobRequest): Operation = withClient { client =>
//    client.deleteJobCallable().call(request)
//  }
//
//  private def withClient[T](f: BatchServiceClient => T): T = {
//    // set user agent to cromwell so requests can be differentiated on batch
//    val headers = ImmutableMap.of("user-agent", "cromwell")
//    val headerProvider = FixedHeaderProvider.create(headers)
//    val batchSettings = BatchServiceSettings.newBuilder.setHeaderProvider(headerProvider).build
//    val client = BatchServiceClient.create(batchSettings)
//    try
//      f(client)
//    finally
//      client.close()
//  }
//}
