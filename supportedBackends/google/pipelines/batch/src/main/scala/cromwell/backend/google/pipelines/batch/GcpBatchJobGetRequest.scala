package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.{BatchServiceClient, JobName}

class GcpBatchJobGetRequest {

  def GetJob(jobName: String) = {

    val projectId = "batch-testing-350715"
    val region = "us-central1"

    val batchServiceClient = BatchServiceClient
      .create()

    var status = "NA"
    println(status)
    while (status != "SUCCEEDED") {
      val job = batchServiceClient
        .getJob(JobName
          .newBuilder()
          .setProject(projectId)
          .setLocation(region)
          .setJob(jobName)
          .build())


      println(job
        .getStatus
        .getState
        .toString)

      status = job
        .getStatus
        .getState
        .toString
      println(f"status in while $status")
    }
    }



}


