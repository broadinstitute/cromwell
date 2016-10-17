package cromwell.backend.impl.tes

import akka.actor.Props

import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import cromwell.backend.wdl.OnlyPureFunctions
import net.ceedubs.ficus.Ficus._

import scalaj.http._

import scala.collection.JavaConverters._
import scala.concurrent.duration.Duration
import scala.concurrent.{Await, Future}

class TesJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                           override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  private val tesEndpoint = configurationDescriptor.backendConfig.as[String]("endpoint")
  private val tesTaskDesc = configurationDescriptor.backendConfig.as[String]("taskDesc")
  private val pollInterval = configurationDescriptor.backendConfig.as[Int]("poll-interval")

  private val tag = s"TesJobExecutionActor-${jobDescriptor.call.fullyQualifiedName}:"
  

  /**
    * Submit a job.
    */
  private def execute(): Future[BackendJobExecutionResponse] = {
    val response: HttpResponse[String] = Http.postData(tesUrl, taskDesc).method("POST")
    val resonseCode = response.code
    log.debug("{} Return code of POST request to TES: {}", tag, responseCode)
    response match {
      case r if r.isSuccess =>
        log.info("{} {} submitted to TES. Waiting for the job to complete.", tag, jobDescriptor.call.fullyQualifiedName)
        val jobId = Json.parse[Map[String,String]](response.body).getOrElse("jobId", None)
        log.debug(s"{} Output of submit process : {}", tag, response.body)
        if (job.nonEmpty) {
          log.info("{} {} mapped to Tes JobID: {}", tag, jobDescriptor.call.fullyQualifiedName, jobId)
          waitUntilDone(jobId)
          SucceededResponse(jobDescriptor.key, Option(0), Map.empty, None, Seq.empty)
        } else {
          self ! JobExecutionResponse(FailedNonRetryableResponse(jobDescriptor.key,
            new IllegalStateException("Failed to retrieve jobId"), Option(response.body)))
        }
      case r if r.isError =>
        self ! JobExecutionResponse(FailedNonRetryableResponse(jobDescriptor.key,
          new IllegalStateException(s"Execution process failed. Tes returned an error: $responseCode"), Option(response.body)))
    }
  }

  private def waitUntilDone(jobId: String): Unit = {
    val response: HttpResponse[String] = Http(endpoint.concat(jobId)).method("GET")

    response match {
      case r if r.isSuccess =>
        val responseMap = Json.parse[Map[String,String]](response.body)
        val statusString = responseMap("state")
        if (statusString == "Complete") {
          log.info(s"Job {} is Complete", jobId)
        } else {
          log.info(s"Still waiting for completion. Last known status: {}", statusString)
          Thread.sleep(pollInterval)
          waitUntilDone(jobId)
        }

      case r if r.isError =>
        val msg = "Could not retreive status from the queue: " + response.body.toString
        logger.error(msg)
        throw new IllegalStateException(msg)
    }
  }
}

object TesJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor,
            configurationDescriptor: BackendConfigurationDescriptor): Props = Props(new TesJobExecutionActor(jobDescriptor, configurationDescriptor))
}
