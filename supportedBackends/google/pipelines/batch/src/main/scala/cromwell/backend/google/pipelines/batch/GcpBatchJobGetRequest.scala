package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.{BatchServiceClient, GetJobRequest, JobName}

class GcpBatchJobGetRequest {

  def GetJob(jobName: String) = {

    val projectId = "batch-testing-350715"
    val region = "us-central1"

    val batchServiceClient = BatchServiceClient.create

    val request = GetJobRequest.newBuilder.setName(JobName.of(projectId, region, jobName).toString()).build
    val job = batchServiceClient.getJob(request)

    batchServiceClient.close()
    job.getStatus.getState
  }

}


