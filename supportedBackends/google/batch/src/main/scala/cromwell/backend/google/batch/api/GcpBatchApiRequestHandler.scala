package cromwell.backend.google.batch.api

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.gax.rpc.FixedHeaderProvider
import com.google.cloud.batch.v1._
import com.google.common.collect.ImmutableMap
import com.google.longrunning.Operation
import common.util.Backoff
import cromwell.backend.google.batch.api.BatchApiRequestManager._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.WorkflowId
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration.DurationInt
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.Try

// A batch of requests to invoke GCP Batch service
//
// GCP Batch does not support batching many requests, we are simulating a way to do it.
// NOTE: This is not thread-safe
// TODO: Alex - what about the performance? PAPIv2 is backed by an ArrayList
//case class GcpBatchGroupedRequests(private val requests: List[BatchApiRequest]) {
//  // prepend operation is more efficient
//  def +(that: BatchApiRequest): GcpBatchGroupedRequests = copy(requests = that :: requests)
//
//  // preserve insertion order
//  def getRequests: List[BatchApiRequest] = requests.reverse
//}
class GcpBatchGroupedRequests {
  private var _requests: List[BatchApiRequest] = List.empty
  def add(that: BatchApiRequest): Unit = _requests = that :: _requests
  def requests: List[BatchApiRequest] = _requests.reverse

  def size: Int = _requests.size

  def execute(): Unit =
    ???
}
/*
public void execute() throws IOException {
  boolean retryAllowed;
  Preconditions.checkState(!requestInfos.isEmpty(), "Batch is empty");

  // Log a warning if the user is using the global batch endpoint. In the future, we can turn this
  // into a preconditions check.
  if (GLOBAL_BATCH_ENDPOINT.equals(this.batchUrl.toString())) {
    LOGGER.log(Level.WARNING, GLOBAL_BATCH_ENDPOINT_WARNING);
  }

  HttpRequest batchRequest = requestFactory.buildPostRequest(this.batchUrl, null);
  // NOTE: batch does not support gzip encoding
  HttpExecuteInterceptor originalInterceptor = batchRequest.getInterceptor();
  batchRequest.setInterceptor(new BatchInterceptor(originalInterceptor));
  int retriesRemaining = batchRequest.getNumberOfRetries();

  do {
    retryAllowed = retriesRemaining > 0;
    MultipartContent batchContent = new MultipartContent();
    batchContent.getMediaType().setSubType("mixed");
    int contentId = 1;
    for (RequestInfo<?, ?> requestInfo : requestInfos) {
      batchContent.addPart(
          new MultipartContent.Part(
              new HttpHeaders().setAcceptEncoding(null).set("Content-ID", contentId++),
              new HttpRequestContent(requestInfo.request)));
    }
    batchRequest.setContent(batchContent);
    HttpResponse response = batchRequest.execute();
    BatchUnparsedResponse batchResponse;
    try {
      // Find the boundary from the Content-Type header.
      String boundary = "--" + response.getMediaType().getParameter("boundary");

      // Parse the content stream.
      InputStream contentStream = new BufferedInputStream(response.getContent());
      batchResponse =
          new BatchUnparsedResponse(contentStream, boundary, requestInfos, retryAllowed);

      while (batchResponse.hasNext) {
        batchResponse.parseNextResponse();
      }
    } finally {
      response.disconnect();
    }

    List<RequestInfo<?, ?>> unsuccessfulRequestInfos = batchResponse.unsuccessfulRequestInfos;
    if (!unsuccessfulRequestInfos.isEmpty()) {
      requestInfos = unsuccessfulRequestInfos;
    } else {
      break;
    }
    retriesRemaining--;
  } while (retryAllowed);
  requestInfos.clear();
}
 */

object GcpBatchGroupedRequests {
  def execute(grouped: GcpBatchGroupedRequests): Unit = {
    grouped.requests.map { request =>
      // this already handles retries
      request.httpRequest.execute()
    }
    ???
  }
}

class GcpBatchApiRequestHandler {
  def submit(request: CreateJobRequest): Job = withClient { client =>
    client.createJobCallable
      .call(request)
  }

  // TODO: Similar to PipelinesApiRequestHandler, this is supposed to not execute the request just yet but enqueue it
  def enqueue(request: BatchApiRequest, batchRequest: GcpBatchGroupedRequests, pollingManager: ActorRef)(implicit
    ec: ExecutionContext
  ): Future[Try[Unit]] = request match {
    case r: BatchRunCreationRequest =>
      val completionPromise = Promise[Try[Unit]]()
      val resultHandler = runCreationResultHandler(r, completionPromise, pollingManager)
      addRunCreationToBatch(r, batchRequest, resultHandler)
      batchRequest.add(r)
      completionPromise.future
  }

  private def runCreationResultHandler(request: BatchRunCreationRequest,
                                       completionPromise: Promise[Try[Unit]],
                                       pollingManager: ActorRef
  ) = {}

  def makeBatchRequest: GcpBatchGroupedRequests = new GcpBatchGroupedRequests

  def query(request: GetJobRequest): Job = withClient { client =>
    client.getJob(request)
  }

  def abort(request: DeleteJobRequest): Operation = withClient { client =>
    client.deleteJobCallable().call(request)
  }

  private def withClient[T](f: BatchServiceClient => T): T = {
    // set user agent to cromwell so requests can be differentiated on batch
    val headers = ImmutableMap.of("user-agent", "cromwell")
    val headerProvider = FixedHeaderProvider.create(headers)
    val batchSettings = BatchServiceSettings.newBuilder.setHeaderProvider(headerProvider).build
    val client = BatchServiceClient.create(batchSettings)
    try
      f(client)
    finally
      client.close()
  }
}
