package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.{BatchServiceClient, JobName}

class GcpBatchStatus {

  val projectId = "batch-testing-350715"
  val region = "us-central1"
  val jobName = "hello-dan-04"

  val batchServiceClient = BatchServiceClient
    .create()

  val job = batchServiceClient
    .getJob(JobName
      .newBuilder()
      .setProject(projectId)
      .setLocation(region)
      .setJob(jobName)
      .build())

}
