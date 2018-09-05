package centaur.api

import java.io.IOException
import java.util.concurrent.Executors

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.model.{HttpRequest, StatusCodes}
import akka.http.scaladsl.unmarshalling.Unmarshaller.UnsupportedContentTypeException
import akka.stream.{ActorMaterializer, ActorMaterializerSettings, BufferOverflowException, StreamTcpException}
import cats.effect.IO
import centaur.test.workflow.Workflow
import centaur.{CentaurConfig, CromwellManager}
import com.typesafe.config.ConfigFactory
import cromwell.api.CromwellClient
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import cromwell.api.model.{CallCacheDiff, CromwellBackends, ShardIndex, SubmittedWorkflow, WorkflowId, WorkflowMetadata, WorkflowOutputs, WorkflowStatus}
import net.ceedubs.ficus.Ficus._

import scala.concurrent._
import scala.concurrent.duration._
import scala.util.Try

object CentaurCromwellClient {
  val LogFailures = ConfigFactory.load().as[Option[Boolean]]("centaur.log-request-failures").getOrElse(false)
  // Do not use scala.concurrent.ExecutionContext.Implicits.global as long as this is using Await.result
  // See https://github.com/akka/akka-http/issues/602
  // And https://github.com/viktorklang/blog/blob/master/Futures-in-Scala-2.12-part-7.md
  final implicit val blockingEc: ExecutionContextExecutor = ExecutionContext.fromExecutor(
    Executors.newFixedThreadPool(100, DaemonizedDefaultThreadFactory))

  // Akka HTTP needs both the actor system and a materializer
  final implicit val system = ActorSystem("centaur-acting-like-a-system")
  final implicit val materializer: ActorMaterializer = ActorMaterializer(ActorMaterializerSettings(system))
  final val apiVersion = "v1"
  val cromwellClient = new CromwellClient(CentaurConfig.cromwellUrl, apiVersion)

  def submit(workflow: Workflow): IO[SubmittedWorkflow] = {
    sendReceiveFutureCompletion(() => cromwellClient.submit(workflow.toWorkflowSubmission(refreshToken = CentaurConfig.optionalToken)))
  }

  def status(workflow: SubmittedWorkflow): IO[WorkflowStatus] = {
    sendReceiveFutureCompletion(() => cromwellClient.status(workflow.id))
  }

  def abort(workflow: SubmittedWorkflow): IO[WorkflowStatus] = {
    sendReceiveFutureCompletion(() => cromwellClient.abort(workflow.id))
  }

  def outputs(workflow: SubmittedWorkflow): IO[WorkflowOutputs] = {
    sendReceiveFutureCompletion(() => cromwellClient.outputs(workflow.id))
  }

  def callCacheDiff(workflowA: SubmittedWorkflow, callA: String, workflowB: SubmittedWorkflow, callB: String): IO[CallCacheDiff] = {
    sendReceiveFutureCompletion(() => cromwellClient.callCacheDiff(workflowA.id, callA, ShardIndex(None), workflowB.id, callB, ShardIndex(None)))
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

  def metadata(workflow: SubmittedWorkflow): IO[WorkflowMetadata] = metadata(workflow.id)

  def metadata(id: WorkflowId): IO[WorkflowMetadata] = {
    sendReceiveFutureCompletion(() => cromwellClient.metadata(id))
  }

  lazy val backends: Try[CromwellBackends] = Try(Await.result(cromwellClient.backends, CromwellManager.timeout * 2))

  def retryRequest[T](x: () => Future[T], timeout: FiniteDuration): IO[T] = {
    // If Cromwell is known not to be ready, delay the request to avoid requests bound to fail 
    val ioDelay = if (!CromwellManager.isReady) IO.sleep(10.seconds) else IO.unit

    ioDelay.flatMap( _ =>
      // Could probably use IO to do the retrying too. For now use a copyport of Retry from cromwell core. Retry 5 times,
      // wait 5 seconds between retries. Timeout the whole thing using the IO timeout.
      IO.fromFuture(IO(Retry.withRetry(x, Option(5), 5.seconds, isTransient = isTransient)).timeout(timeout))
    )
  }

  def sendReceiveFutureCompletion[T](x: () => Future[T]) = {
    retryRequest(x, CentaurConfig.sendReceiveTimeout)
  }

  private def isTransient(f: Throwable) = f match {
    case _: TimeoutException |
                    _: StreamTcpException |
                    _: IOException |
                    _: UnsupportedContentTypeException => true
    case BufferOverflowException(message) => message.contains("Please retry the request later.")
    case unsuccessful: UnsuccessfulRequestException => unsuccessful.httpResponse.status == StatusCodes.NotFound
    case unexpected: RuntimeException => unexpected.getMessage.contains("The http server closed the connection unexpectedly") 
    case _ => false
  }
}
