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

    println(job.getStatus.getState.toString)

    val status = job.getStatus.getState.toString
    //Future.successful(GcpBatchRunStatus.Running)

    if (status == "SUCCEEDED") {

      //val testStatus = GcpBatchJobSuccess(jobName = jobName, result = status)
      //testStatus
      status
      //Future.successful(GcpBatchRunStatus.Complete)

    }
    else println(f"status is $status")
  }


}


