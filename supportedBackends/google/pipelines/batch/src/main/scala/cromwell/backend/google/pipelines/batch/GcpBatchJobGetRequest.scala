package cromwell.backend.google.pipelines.batch

//import cats.instances.unit
//import com.google.api.core.ApiFuture
//import cromwell.backend.google.pipelines.batch.GcpBatchBackendSingletonActor.BatchGetJob

//import scala.concurrent.Promise
import scala.util.Try
import com.google.cloud.batch.v1.{BatchServiceClient, GetJobRequest, JobName}
//import com.google.cloud.batch.v1.{BatchServiceClient, GetJobRequest, Job, JobName}

//import scala.concurrent.Await
//import scala.concurrent.duration.DurationInt
//import scala.util.{Failure, Success, Try}
//import org.http4s.Response.timeout

//import scala.concurrent.Await
//import scala.concurrent.duration.Duration
//import cromwell.backend.google.pipelines.batch.GcpBatchBackendSingletonActor.GcpBatchJobSuccess

//import scala.concurrent.Future

class GcpBatchJobGetRequest {

  def GetJob(jobName: String) = {

    val projectId = "batch-testing-350715"
    val region = "us-central1"

    val batchServiceClient = BatchServiceClient.create


    val request = GetJobRequest.newBuilder.setName(JobName.of(projectId, region, jobName).toString()).build
    val job = batchServiceClient.getJob(request)

    batchServiceClient.close()
    job





    //print(response.getStatus.getState)
    //print(response.getStatus.toString)
    //batchServiceClient.close()

    //response


    /*
        val batchServiceClient = BatchServiceClient.create()

        val job = batchServiceClient
          .getJob(JobName
            .newBuilder()
            .setProject(projectId)
            .setLocation(region)
            .setJob(jobName)
            .build())
        */


    /*
    val batchServiceClient = BatchServiceClient.create
    val request = GetJobRequest.newBuilder.setName(JobName.of(projectId, region, jobName).toString())
                               .build
    println(request.toString)
    val future: ApiFuture[Job] = batchServiceClient.getJobCallable.futureCall(request)
    //Await.ready(future.get, Duration.Inf)
    //future.get

    val response = future.get()
    //println(response)
    //response.getStatus.getState
    batchServiceClient.close()
    //println(response.getStatus.getState)
    */
    //BatchGetJob(jobName)

    /*
    Await.result(future, 1.second) match {
      case Success(res) => println("succeeded")
      case Failure(e) => println("failed")
    }*/


    //val response: Job = future.get()


    //val status = job.getStatus.getState
    //val jobResult = job
    //jobResult

  }

  def jobGetRequest(jobId: String) = {
    val gcpBatchPoll = new GcpBatchJobGetRequest
    gcpBatchPoll.GetJob(jobId)
    //jobDetail
  }

  def status(jobId: String): Try[RunStatus] = for {
    _ <- Try(jobGetRequest(jobId).toString)
    //runStatus <- RunStatus.fromJobStatus(jobId)
    runStatus <- RunStatus.testJobStatus(jobId)
  } yield runStatus


}


