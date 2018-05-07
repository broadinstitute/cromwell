package centaur.api

import java.io.IOException
import java.util.concurrent.Executors

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.model.{HttpRequest, StatusCodes}
import akka.http.scaladsl.unmarshalling.Unmarshaller.UnsupportedContentTypeException
import akka.stream.{ActorMaterializer, ActorMaterializerSettings, StreamTcpException}
import centaur.test.metadata.WorkflowMetadata
import centaur.test.workflow.Workflow
import centaur.{CentaurConfig, CromwellManager}
import cromwell.api.CromwellClient
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import cromwell.api.model.{CromwellBackends, SubmittedWorkflow, WorkflowId, WorkflowOutputs, WorkflowStatus}

import scala.concurrent._
import scala.concurrent.duration._
import scala.util.Try

object CentaurCromwellClient {
  // Do not use scala.concurrent.ExecutionContext.Implicits.global as long as this is using Await.result
  // See https://github.com/akka/akka-http/issues/602
  // And https://github.com/viktorklang/blog/blob/master/Futures-in-Scala-2.12-part-7.md
  final implicit val blockingEc: ExecutionContextExecutor = ExecutionContext.fromExecutor(
    Executors.newCachedThreadPool(DaemonizedDefaultThreadFactory))

  // Akka HTTP needs both the actor system and a materializer
  final implicit val system = ActorSystem("centaur-acting-like-a-system")
  final implicit val materializer: ActorMaterializer = ActorMaterializer(ActorMaterializerSettings(system))
  final val apiVersion = "v1"
  val cromwellClient = new CromwellClient(CentaurConfig.cromwellUrl, apiVersion)

  def submit(workflow: Workflow): Try[SubmittedWorkflow] = {
    sendReceiveFutureCompletion(() => cromwellClient.submit(workflow.toWorkflowSubmission(refreshToken = CentaurConfig.optionalToken)))
  }

  def status(workflow: SubmittedWorkflow): Try[WorkflowStatus] = {
    sendReceiveFutureCompletion(() => cromwellClient.status(workflow.id))
  }

  def abort(workflow: SubmittedWorkflow): Try[WorkflowStatus] = {
    sendReceiveFutureCompletion(() => cromwellClient.abort(workflow.id))
  }

  def outputs(workflow: SubmittedWorkflow): Try[WorkflowOutputs] = {
    sendReceiveFutureCompletion(() => cromwellClient.outputs(workflow.id))
  }

  /*
    Sends a quick ping to the Cromwell query endpoint. The query endpoint is the only one which both hits the
    database w/o requiring a workflow id and does not modify server state. Not using CromwellClient here as it
    currently does not support query.
   */
  def isAlive: Boolean = {
    val request = Http().singleRequest(HttpRequest(uri=s"${CentaurConfig.cromwellUrl}/api/workflows/$apiVersion/query?status=Succeeded"))
    Try(Await.result(request, CentaurConfig.sendReceiveTimeout)).isSuccess
  }

  def metadata(workflow: SubmittedWorkflow): Try[WorkflowMetadata] = metadata(workflow.id)

  def metadata(id: WorkflowId): Try[WorkflowMetadata] = {
    sendReceiveFutureCompletion(() => cromwellClient.metadata(id)) map { m =>
      WorkflowMetadata.fromMetadataJson(m.value).toOption.get
    }
  }

  lazy val backends: Try[CromwellBackends] = sendReceiveFutureCompletion(() => cromwellClient.backends)

  // The total time waiting for a Future, including network hiccups, must be longer than the time for a restart.
  private val awaitTotalTimeout: FiniteDuration = CromwellManager.timeout * 2
  private val awaitSleep: FiniteDuration = 5.seconds
  private val awaitMaxAttempts: Long = awaitTotalTimeout.toSeconds / awaitSleep.toSeconds

  /**
    * Ensure that the Future completes within the specified timeout. If it does not, or if the Future fails,
    * will return a Failure, otherwise a Success
    */
  def awaitFutureCompletion[T](x: () => Future[T], timeout: FiniteDuration, attempt: Int = 1): Try[T] = {
    // We can't recover the future itself with a "recoverWith retry pattern" because it'll timeout anyway from the Await.result
    // We want to keep timing out to catch cases where Cromwell becomes unresponsive
    Try(Await.result(x(), timeout)) recoverWith {
      case _: TimeoutException |
           _: StreamTcpException |
           _: IOException |
           _: UnsupportedContentTypeException if attempt < awaitMaxAttempts =>
        Thread.sleep(awaitSleep.toMillis)
        awaitFutureCompletion(x, timeout, attempt + 1)
      case unsuccessful: UnsuccessfulRequestException if unsuccessful.httpResponse.status == StatusCodes.NotFound && attempt < awaitMaxAttempts =>
        Thread.sleep(awaitSleep.toMillis)
        awaitFutureCompletion(x, timeout, attempt + 1)
      // see https://github.com/akka/akka-http/issues/768
      case unexpected: RuntimeException
        if unexpected.getMessage.contains("The http server closed the connection unexpectedly") &&
          attempt < awaitMaxAttempts =>
        Thread.sleep(awaitSleep.toMillis)
        awaitFutureCompletion(x, timeout, attempt + 1)
    }
  }

  def sendReceiveFutureCompletion[T](x: () => Future[T]) = {
    awaitFutureCompletion(x, CentaurConfig.sendReceiveTimeout)
  }

  def maxWorkflowLengthCompletion[T](x: () => Future[T]) = {
    awaitFutureCompletion(x, CentaurConfig.maxWorkflowLength)
  }
}
