package centaur.api

import java.io.IOException
import java.util.concurrent.Executors

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.model.StatusCodes.ClientError
import akka.http.scaladsl.model.{HttpRequest, StatusCodes}
import akka.http.scaladsl.unmarshalling.Unmarshaller.UnsupportedContentTypeException
import akka.stream.{ActorMaterializer, ActorMaterializerSettings, BufferOverflowException, StreamTcpException}
import cats.effect.IO
import centaur.test.workflow.Workflow
import centaur.{CentaurConfig, CromwellManager}
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.StrictLogging
import cromwell.api.CromwellClient
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import cromwell.api.model._
import net.ceedubs.ficus.Ficus._

import scala.concurrent._
import scala.concurrent.duration._
import scala.util.Try

object CentaurCromwellClient extends StrictLogging {
  val config = ConfigFactory.load()
  val LogFailures = config.as[Option[Boolean]]("centaur.log-request-failures").getOrElse(false)
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
  
  val defaultMetadataArgs = config
    .getAs[Map[String, List[String]]]("centaur.metadata-args")

  def submit(workflow: Workflow): IO[SubmittedWorkflow] = {
    sendReceiveFutureCompletion(() => {
      val submitted = cromwellClient.submit(workflow.toWorkflowSubmission(refreshToken = CentaurConfig.optionalToken))
      submitted.biSemiflatMap(
        httpResponse =>
          for {
            _ <- httpResponse.status match {
              case _: ClientError => IO(logger.info(s"Submitting ${workflow.testName} returned ${httpResponse.status}"))
              case _ => IO(logger.error(s"Submitting ${workflow.testName} returned unexpected ${httpResponse.status}"))
            }
          } yield httpResponse,
        submittedWorkflow =>
          for {
            _ <- IO(logger.info(s"Submitting ${workflow.testName} returned workflow id ${submittedWorkflow.id}"))
          } yield submittedWorkflow
      )
    })
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
    val response = Http().singleRequest(HttpRequest(uri=s"${CentaurConfig.cromwellUrl}/api/workflows/$apiVersion/query?status=Succeeded"))
    // Silence the following warning by discarding the result of a successful query:
    // Response entity was not subscribed after 1 second. Make sure to read the response entity body or call `discardBytes()` on it.
    val successOrFailure = response map { _.entity.discardBytes() }
    Try(Await.result(successOrFailure, CentaurConfig.sendReceiveTimeout)).isSuccess
  }

  def metadata(workflow: SubmittedWorkflow, args: Option[Map[String, List[String]]] = defaultMetadataArgs): IO[WorkflowMetadata] = metadataWithId(workflow.id, args)

  def metadataWithId(id: WorkflowId, args: Option[Map[String, List[String]]] = defaultMetadataArgs): IO[WorkflowMetadata] = {
    sendReceiveFutureCompletion(() => cromwellClient.metadata(id, args))
  }
  
  implicit private val timer = IO.timer(blockingEc)
  implicit private val contextShift = IO.contextShift(blockingEc)

  lazy val backends: IO[CromwellBackends] = cromwellClient.backends.timeout(CromwellManager.timeout * 2).asIo

  def retryRequest[T](func: () => FailureResponseOrT[T], timeout: FiniteDuration): IO[T] = {
    // If Cromwell is known not to be ready, delay the request to avoid requests bound to fail
    val ioDelay = if (!CromwellManager.isReady) IO.sleep(10.seconds) else IO.unit

    ioDelay.flatMap( _ =>
      // Could probably use IO to do the retrying too. For now use a copyport of Retry from cromwell core. Retry 5 times,
      // wait 5 seconds between retries. Timeout the whole thing using the IO timeout.
      // https://github.com/cb372/cats-retry
      // https://typelevel.org/cats-effect/datatypes/io.html#example-retrying-with-exponential-backoff
      IO.fromFuture(IO(Retry.withRetry(
        () => func().asIo.unsafeToFuture(),
        Option(5),
        5.seconds,
        isTransient = isTransient)
      ).timeout(timeout))
    )
  }

  def sendReceiveFutureCompletion[T](x: () => FailureResponseOrT[T]): IO[T] = {
    retryRequest(x, CentaurConfig.sendReceiveTimeout)
  }

  private def isTransient(f: Throwable) = f match {
    case _: StreamTcpException |
         _: IOException |
         _: UnsupportedContentTypeException => true
    case BufferOverflowException(message) => message.contains("Please retry the request later.")
    case unsuccessful: UnsuccessfulRequestException => unsuccessful.httpResponse.status == StatusCodes.NotFound
    case unexpected: RuntimeException => unexpected.getMessage.contains("The http server closed the connection unexpectedly") 
    case _ => false
  }
}
