package cromwell.backend.impl.tes

import akka.actor.{Actor, ActorLogging, ActorSystem, Props}
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.async.AsyncBackendJobExecutionActor.ExecutionMode
import cromwell.backend.async.{AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, NonRetryableExecution, SuccessfulExecutionHandle}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.backend.impl.tes.util._
import cromwell.core.logging.JobLogging
import cromwell.core.retry.SimpleExponentialBackoff
import TesResponseJsonFormatter._
import spray.httpx.SprayJsonSupport._
import spray.client.pipelining._
import spray.http.HttpRequest
import spray.httpx.unmarshalling._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.language.postfixOps

final case class TesAsyncBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                                  override val completionPromise: Promise[BackendJobExecutionResponse],
                                                  configurationDescriptor: BackendConfigurationDescriptor)
  extends Actor with ActorLogging with AsyncBackendJobExecutionActor with JobLogging {
  import TesAsyncBackendJobExecutionActor._

  private implicit val actorSystem: ActorSystem = context.system
  private def pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]

  override lazy val jobTag = jobDescriptor.key.tag
  override lazy val workflowId = jobDescriptor.workflowDescriptor.id
  override lazy val pollBackOff = SimpleExponentialBackoff(initialInterval = 1 seconds, maxInterval = 10 seconds, multiplier = 1.1)
  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(initialInterval = 1 seconds, maxInterval = 20 seconds, multiplier = 1.1)
  override lazy val retryable = false

  private val tesEndpoint = configurationDescriptor.backendConfig.as[String]("endpoint")

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    previous match {
      case handle: TesPendingExecutionHandle =>
        jobLogger.info(s"Polling TES Job ${handle.job.jobId}")
        updateExecutionHandle(handle)
      case f: FailedNonRetryableExecutionHandle => f.future
      case s: SuccessfulExecutionHandle => s.future
      case badHandle => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $badHandle"))
    }
  }

  private def updateExecutionHandle(oldHandle: TesPendingExecutionHandle): Future[ExecutionHandle] = {
    def successulResponse(response: TesGetResponse): ExecutionHandle = {
      if (response.state contains "Complete") {
        jobLogger.info(s"Job ${oldHandle.job.jobId} is complete")
        SuccessfulExecutionHandle(Map.empty, 0, Map.empty, Seq.empty, None) // FIXME: blah
      } else {
        jobLogger.info(s"Status is ${response.state}")
        oldHandle
      }
    }

    pipeline[TesGetResponse].apply(Get(s"$tesEndpoint/${oldHandle.job.jobId}")) map successulResponse recover failedTesResponse
  }

  def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    def successfulResponse(response: TesPostResponse): ExecutionHandle = {
      response.value match {
        case Some(jobId) =>
          jobLogger.info(s"Launched $workflowId as TES job $jobId")
          TesPendingExecutionHandle(jobDescriptor, TesJob(jobId))
        case None =>
          FailedNonRetryableExecutionHandle(new Exception(s"Unable to retrieve TES job ID for $workflowId"))
      }
    }

    // FIXME: Only executing now, no recover
    pipeline[TesPostResponse].apply(Post(tesEndpoint, TaskDesc)) map successfulResponse recover failedTesResponse
  }

  override protected implicit def ec: ExecutionContext = context.dispatcher
}

object TesAsyncBackendJobExecutionActor {
  def failedTesResponse: PartialFunction[Throwable, ExecutionHandle] = {
    case e => FailedNonRetryableExecutionHandle(e)
  }

  def props(jobDescriptor: BackendJobDescriptor,
            completionPromise: Promise[BackendJobExecutionResponse],
            configurationDescriptor: BackendConfigurationDescriptor): Props = {
    Props(TesAsyncBackendJobExecutionActor(jobDescriptor, completionPromise, configurationDescriptor))
  }

  final case class TesPendingExecutionHandle(jobDescriptor: BackendJobDescriptor, job: TesJob) extends ExecutionHandle {
    override val isDone = false
    override val result = NonRetryableExecution(new IllegalStateException("TesPendingExecutionHandle cannot yield a result"))
  }

  final case class TesJob(jobId: String)

  val TaskDesc = TesTask( // FIXME
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