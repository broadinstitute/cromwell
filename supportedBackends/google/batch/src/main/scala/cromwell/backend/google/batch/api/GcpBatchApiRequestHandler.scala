package cromwell.backend.google.batch.api

import com.google.api.gax.rpc.FixedHeaderProvider
import com.google.cloud.batch.v1._
import com.google.common.collect.ImmutableMap
import com.google.longrunning.Operation

class GcpBatchApiRequestHandler {
  def submit(request: CreateJobRequest): Job = withClient { client =>
    client.createJobCallable
      .call(request)
  }

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
