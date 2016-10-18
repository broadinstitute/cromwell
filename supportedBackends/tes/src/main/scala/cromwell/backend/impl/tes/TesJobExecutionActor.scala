package cromwell.backend.impl.tes

import akka.actor.Props

import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse, FailedNonRetryableResponse}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import net.ceedubs.ficus.Ficus._
import cromwell.backend.impl.tes.util._
import TesClient._
import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.util.{Success, Failure, Try}

class TesJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                           override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  val tesEndpoint = configurationDescriptor.backendConfig.as[String]("endpoint")
  val tesTaskDesc = configurationDescriptor.backendConfig.as[String]("taskDesc")
  val pollInterval = configurationDescriptor.backendConfig.as[Long]("poll-interval")

  private val tag = s"TesJobExecutionActor-${jobDescriptor.call.fullyQualifiedName}:"
 
  /**
    * Submit a job.
    */
  override def execute: Future[BackendJobExecutionResponse] = {
    val response = Try(Await.result(Pipeline[TesPostResponse].apply(Post(tesEndpoint, tesTaskDesc)), 5.seconds))

    response match {
      case Success(r) =>
        log.info("{} {} submitted to TES. Waiting for the job to complete.", tag, jobDescriptor.call.fullyQualifiedName)
        val jobId: String = r.value.get // FIXME: remove the .get
        if (jobId.nonEmpty) {
          log.info("{} {} mapped to Tes JobID: {}", tag, jobDescriptor.call.fullyQualifiedName, jobId)
          waitUntilDone(jobId)
          Future.successful(SucceededResponse(jobDescriptor.key, Option(0), Map.empty, None, Seq.empty))
        } else {
          Future.successful(FailedNonRetryableResponse(jobDescriptor.key,
            new IllegalStateException("Failed to retrieve jobId"), Option(1)))
        }

      case Failure(e) =>
        Future.successful(FailedNonRetryableResponse(jobDescriptor.key,
          new IllegalStateException(s"Execution process failed. Tes returned an error: ${e.getMessage}"), Option(1)))
    }
  }

  private def waitUntilDone(jobId: String): Unit = {
    val response = Try(Await.result(Pipeline[TesGetResponse].apply(Get(s"$tesEndpoint/$jobId")), 5.seconds))
    
    response match {
      case Success(r) =>
        val statusString = r.state
        if (statusString.contains("Complete")) {
          log.info(s"Job {} is Complete", jobId)
        } else {
          log.info(s"Still waiting for completion. Last known status: {}", statusString)
          Thread.sleep(pollInterval)
          waitUntilDone(jobId)
        }

      case Failure(e) =>
        val msg = s"Could not retreive status from the queue: ${e.getMessage}"
        throw new IllegalStateException(msg)
    }
  }
}

object TesJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor,
            configurationDescriptor: BackendConfigurationDescriptor): Props = Props(new TesJobExecutionActor(jobDescriptor, configurationDescriptor))

  val foo = TesTask(
    Some("TestMD5"),
    Some("MyProject"),
    Some("My Desc"),
    None,
    Some(Seq(TaskParameter(None, None, None, Some("/tmp/test_out"), None, None))),
    Some(Resources(None, None, None, Some(Seq(Volume(Some("test_file"), Some(1), None, Some("/tmp")))), None)),
    None,
    Some(Seq(DockerExecutor(Some("ubuntu"), Some(Seq("echo", "foo")), None, Some("/tmp/test_out"), None)))
  )
}
