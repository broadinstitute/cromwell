package centaur.test

import java.util.UUID

import akka.actor.ActorSystem
import centaur.api.CromwellStatusJsonSupport
import cats.Monad
import centaur._
import centaur.api.{CromwellStatus, _}
import spray.client.pipelining._
import spray.http.{FormData, HttpRequest, HttpResponse}
import spray.httpx.PipelineException
import spray.httpx.unmarshalling._
import spray.json.{DefaultJsonProtocol, JsArray, JsString, JsValue}

import scala.annotation.tailrec
import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{Await, Future}
import scala.util.{Failure, Success, Try}
import spray.httpx.SprayJsonSupport._
import FailedWorkflowSubmissionJsonSupport._
import CromwellStatusJsonSupport._
import centaur.test.metadata.WorkflowMetadata
import centaur.test.workflow.Workflow

import scala.concurrent.ExecutionContext.Implicits.global

/**
  * A simplified riff on the final tagless pattern where the interpreter (monad & related bits) are fixed. Operation
  * functions create an instance of a Test and override the run method to do their bidding. It is unlikely that you
  * should be modifying Test directly, instead most likely what you're looking to do is add a function to the Operations
  * object below
  */
sealed abstract class Test[A] {
  def run: Try[A]
}

object Test {
  implicit val testMonad: Monad[Test] = new Monad[Test] {
    override def flatMap[A, B](fa: Test[A])(f: A => Test[B]): Test[B] = {
      new Test[B] {
        override def run: Try[B] = fa.run flatMap { f(_).run }
      }
    }

    override def pure[A](x: A): Test[A] = {
      new Test[A] {
        override def run: Try[A] = Try(x)
      }
    }
  }
}

/**
  * Defines functions which are building blocks for test formulas. Each building block is expected to perform
  * a single task and these tasks can be composed together to form arbitrarily complex test strategies. For instance
  * submitting a workflow, polling until a status is reached, retrieving metadata, verifying some value, delaying for
  * N seconds - these would all be operations.
  *
  * All operations are expected to return a Test type and implement the run method. These can then
  * be composed together via a for comprehension as a test formula and then run by some other entity.
  */
object Operations {
  def submitWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = {
    new Test[SubmittedWorkflow] {
      override def run: Try[SubmittedWorkflow] = {
        // Collect only the parameters which exist:
        val params = List("wdlSource" -> Option(workflow.data.wdl),
          "workflowInputs" -> workflow.data.inputs,
          "workflowOptions" -> workflow.data.options) collect { case (name, Some(value)) => (name, value) }
        val formData = FormData(params)
        val response = Pipeline[CromwellStatus].apply(Post(CentaurConfig.cromwellUrl + "/api/workflows/v1", formData))
        sendReceiveFutureCompletion(response map { _.id } map UUID.fromString map { SubmittedWorkflow(_, CentaurConfig.cromwellUrl, workflow) })
      }
    }
  }

  /**
    * Polls until a specific status is reached. If a terminal status which wasn't expected is returned, the polling
    * stops with a failure.
    */
  def pollUntilStatus(workflow: SubmittedWorkflow, expectedStatus: WorkflowStatus): Test[SubmittedWorkflow] = {
    new Test[SubmittedWorkflow] {
      @tailrec
      def doPerform(): SubmittedWorkflow = {
        val response = Pipeline[CromwellStatus].apply(Get(CentaurConfig.cromwellUrl + "/api/workflows/v1/" + workflow.id + "/status"))
        val status = sendReceiveFutureCompletion(response map { r => WorkflowStatus(r.status) })
        status match {
          case Success(s) if s == expectedStatus => workflow
          case Success(s: TerminalStatus) => throw new Exception(s"Unexpected terminal status $s but was waiting for $expectedStatus")
          case Failure(f) => throw f
          case _ =>
            Thread.sleep(10000) // This could be a lot smarter including cromwell style backoff
            doPerform()
        }
      }

      override def run: Try[SubmittedWorkflow] = workflowLengthFutureCompletion(Future { doPerform() })
    }
  }

  def validateMetadata(workflow: SubmittedWorkflow, expectedMetadata: WorkflowMetadata): Test[Unit] = {
    val consistencyTimeout = Future { Thread.sleep(CentaurConfig.metadataConsistencyTimeout.toMillis) }

    new Test[Unit] {
      def validateMetadataUntilTimeout(workflow: SubmittedWorkflow, expectedMetadata: WorkflowMetadata): Try[Unit] = {
        def checkDiff(diffs: Iterable[String]): Try[Unit] = {
          diffs match {
            case d if d.isEmpty => Success()
            case d if consistencyTimeout.isCompleted =>
              Failure(throw new Exception(s"Invalid metadata response:\n -${d.mkString("\n -")}\n"))
            case _ => validateMetadataUntilTimeout(workflow, expectedMetadata)
          }
        }

        val response = MetadataRequest(Get(CentaurConfig.cromwellUrl + "/api/workflows/v1/" + workflow.id + "/metadata"))
        // Try to convert the response to a Metadata in our return Try.
        // Currently any error msg will be opaque as it's unlikely to be an issue (can flesh out later)
        val metadata = sendReceiveFutureCompletion(response map { r => WorkflowMetadata.fromMetadataJson(r).toOption.get })
        metadata.map(expectedMetadata.diff(_, workflow.id)) flatMap checkDiff
      }

      override def run: Try[Unit] = {
        validateMetadataUntilTimeout(workflow, expectedMetadata)
      }
    }
  }

  /**
    * Verify that none of the calls within the workflow are cached.
    */
  def validateCachingWasOff(metadata: WorkflowMetadata, workflowName: String): Test[Unit] = {
    new Test[Unit] {
      override def run: Try[Unit] = {
        val cacheHits = metadata.value.keySet filter { _.contains("cacheHitCall") }

        if (cacheHits.isEmpty) Success(())
        else Failure(new Exception(s"Found unexpected cache hits for $workflowName: \n"))
      }
    }
  }

  // FIXME: Should be abstracted w/ validateMetadata - ATM still used by the unused caching tests
  def retrieveMetadata(workflow: SubmittedWorkflow): Test[WorkflowMetadata] = {
    new Test[WorkflowMetadata] {
      override def run: Try[WorkflowMetadata] = {
        val response = MetadataRequest(Get(CentaurConfig.cromwellUrl + "/api/workflows/v1/" + workflow.id + "/metadata"))
        // Try to convert the response to a Metadata in our return Try.
        // Currently any error msg will be opaque as it's unlikely to be an issue (can flesh out later)
        sendReceiveFutureCompletion(response map { r => WorkflowMetadata.fromMetadataJson(r).toOption.get })
      }
    }
  }

  /**
    * Ensure that the Future completes within the specified timeout. If it does not, or if the Future fails,
    * will return a Failure, otherwise a Success
    */
  def awaitFutureCompletion[T](x: Future[T], timeout: FiniteDuration) = Try(Await.result(x, timeout))
  def sendReceiveFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.sendReceiveTimeout)
  def workflowLengthFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.maxWorkflowLength)
  def metadataFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.metadataConsistencyTimeout)

  // Spray needs an implicit ActorSystem
  implicit val system = ActorSystem("centaur-foo")
  val MetadataRequest: HttpRequest => Future[String] = sendReceive ~> unmarshal[String]

  def Pipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[T] = sendReceive ~> unmarshal[T]
  def FailingPipeline[T: FromResponseUnmarshaller]: HttpRequest => Future[FailedWorkflowSubmission] = sendReceive ~> unmarshalFailure[FailedWorkflowSubmission]

  private def unmarshalFailure[T: FromResponseUnmarshaller]: HttpResponse => T =
    response =>
      if (!response.status.isSuccess)
        response.as[T] match {
          case Right(value) => value
          case Left(error: MalformedContent) =>
            throw new PipelineException(error.errorMessage, error.cause.orNull)
          case Left(error) => throw new PipelineException(error.toString)
        }
      else throw new SuccessfulResponseException(response)

  class SuccessfulResponseException(val response: HttpResponse) extends RuntimeException(s"Status: ${response.status}\n" +
    s"Body: ${if (response.entity.data.length < 1024) response.entity.asString else response.entity.data.length + " bytes"}")
}
