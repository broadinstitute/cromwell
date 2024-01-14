package cromwell.backend.google.batch.api

import akka.actor.ActorRef
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.gax.rpc.FixedHeaderProvider
import com.google.cloud.batch.v1._
import com.google.common.collect.ImmutableMap
import com.google.longrunning.Operation

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

// TODO: Similar to PAPIApiRequest
// Lets add the supported request types
sealed trait BatchApiRequest

class GcpBatchApiRequestHandler {
  def submit(request: CreateJobRequest): Job = withClient { client =>
    client.createJobCallable
      .call(request)
  }

  // TODO: Similar to PipelinesApiRequestHandler, this is supposed to not execute the request just yet but enqueue it
  def enqueue(request: BatchApiRequest, batchRequest: BatchRequest, pollingManager: ActorRef)(implicit
    ec: ExecutionContext
  ): Future[Try[Unit]] = ???

  // TODO: Similar to PipelinesApiRequestHandler, this should create a request where we can add many operations to the batch
  // it is likely that we need to use a different model because BatchRequest is deprecated
  def makeBatchRequest: BatchRequest = ???

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
