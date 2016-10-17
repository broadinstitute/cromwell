package cromwell.backend.impl.tes

import akka.actor.Props

import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse, FailedNonRetryableResponse}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import net.ceedubs.ficus.Ficus._
import spray.json._
import cromwell.backend.impl.tes.util._
import TesResponseJsonFormatter._
import TesClient._
import scala.concurrent.duration._
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.{Await, Future}
import scala.util.{Success, Failure, Try}

class TesJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                           override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  private val tesEndpoint = configurationDescriptor.backendConfig.as[String]("endpoint")
  private val tesTaskDesc = configurationDescriptor.backendConfig.as[String]("taskDesc")
  private val pollInterval = configurationDescriptor.backendConfig.as[Long]("poll-interval")

  private val tag = s"TesJobExecutionActor-${jobDescriptor.call.fullyQualifiedName}:"
 
  /**
    * Submit a job.
    */
  override def execute: Future[BackendJobExecutionResponse] = {
    val response = Try(Await.result(Pipeline[TesGetResponse].apply(Post(tesEndpoint, tesTaskDesc)), 5.seconds))

    response match {
      case Success =>
        log.info("{} {} submitted to TES. Waiting for the job to complete.", tag, jobDescriptor.call.fullyQualifiedName)
        val jobId: String = response.value.get // FIXME: remove the .get
        log.debug(s"{} Output of submit process : {}", tag, response.body)
        if (jobId.nonEmpty) {
          log.info("{} {} mapped to Tes JobID: {}", tag, jobDescriptor.call.fullyQualifiedName, jobId)
          waitUntilDone(jobId)
          Future.successful(SucceededResponse(jobDescriptor.key, Option(0), Map.empty, None, Seq.empty))
        } else {
          Future.successful(FailedNonRetryableResponse(jobDescriptor.key,
            new IllegalStateException("Failed to retrieve jobId"), Option(1)))
        }

      case Failure =>
        Future.successful(FailedNonRetryableResponse(jobDescriptor.key,
          new IllegalStateException(s"Execution process failed. Tes returned an error."), Option(1)))
    }
  }

  private def waitUntilDone(jobId: String): Unit = {
    val response = Try(Await.result(Pipeline[TesGetResponse].apply(Get(s"$tesEndpoint/$jobId")), 5.seconds))
    
    response match {
      case Success =>
        val statusString = response.state
        if (statusString.contains("Complete")) {
          log.info(s"Job {} is Complete", jobId)
        } else {
          log.info(s"Still waiting for completion. Last known status: {}", statusString)
          Thread.sleep(pollInterval)
          waitUntilDone(jobId)
        }

      case Failure =>
        val msg = "Could not retreive status from the queue."
        throw new IllegalStateException(msg)
    }
  }
}

object TesJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor,
            configurationDescriptor: BackendConfigurationDescriptor): Props = Props(new TesJobExecutionActor(jobDescriptor, configurationDescriptor))
}
