package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.{BatchServiceClient, JobName}
//import cromwell.backend.google.pipelines.batch.GcpBatchBackendSingletonActor.GcpBatchJobSuccess

//import scala.concurrent.Future

final class GcpBatchJobGetRequest {

  def GetJob(jobName: String) = {

    val projectId = "batch-testing-350715"
    val region = "us-central1"

    val batchServiceClient = BatchServiceClient.create()

    val job = batchServiceClient
      .getJob(JobName
        .newBuilder()
        .setProject(projectId)
        .setLocation(region)
        .setJob(jobName)
        .build())

    Thread.sleep(10000) // testing to see how it affects job polling and if job shows up

    //val status = job.getStatus.getState
    val jobResult = job
    jobResult

  }


}


