package centaur

import java.util.UUID
import java.nio.file.Files

import akka.actor.ActorSystem
import centaur.api._
import centaur.api.{CromwellStatusJsonSupport, CromwellStatus}
import cats.Monad
import spray.client.pipelining._
import spray.http.{HttpResponse, HttpRequest, FormData}
import spray.httpx.{PipelineException, UnsuccessfulResponseException}
import spray.httpx.unmarshalling._
import better.files._

import scala.annotation.tailrec
import scala.concurrent.{Await, Future}
import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success, Try}
import spray.httpx.SprayJsonSupport._
import FailedWorkflowSubmissionJsonSupport._
import CromwellStatusJsonSupport._
import spray.json._

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
  def submitWorkflow(request: WorkflowRequest): Test[Workflow] = {
    new Test[Workflow] {
      override def run: Try[Workflow] = {
        // Collect only the parameters which exist:
        val params = List("wdlSource" -> Option(request.wdl), "workflowInputs" -> request.inputs, "workflowOptions" -> request.options) collect {
          case (name, Some(value)) => (name, value)
        }
        val formData = FormData(params)
        val response = Pipeline[CromwellStatus].apply(Post(CentaurConfig.cromwellUrl + "/api/workflows/v1", formData))
        sendReceiveFutureCompletion(response map { _.id } map UUID.fromString map { Workflow(_, CentaurConfig.cromwellUrl) })
      }
    }
  }

  /**
    * Submits a workflow and expects the response from the server to be an error code.
    *
    * @return The message in the error response
    */
  def submitWorkflowExpectingRejection(request: WorkflowRequest): Test[String] = {
    new Test[String] {
      override def run: Try[String] = {
        // Collect only the parameters which exist:
        val params = List("wdlSource" -> Option(request.wdl), "workflowInputs" -> request.inputs, "workflowOptions" -> request.options) collect {
          case (name, Some(value)) => (name, value)
        }
        val formData = FormData(params)
        val response = FailingPipeline[FailedWorkflowSubmission].apply(Post(CentaurConfig.cromwellUrl + "/api/workflows/v1", formData))
        sendReceiveFutureCompletion(response map { _.message })
      }
    }
  }

  /**
    * Polls until a specific status is reached. If a terminal status which wasn't expected is returned, the polling
    * stops with a failure.
    */
  def pollUntilStatus(workflow: Workflow, expectedStatus: WorkflowStatus): Test[Workflow] = {
    new Test[Workflow] {
      @tailrec
      def doPerform(): Workflow = {
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

      override def run: Try[Workflow] = workflowLengthFutureCompletion(Future { doPerform() } )
    }
  }

  /**
    * This method uses a Metadata JSON returned by Cromwell as it's parameter,
    * and creates a test map. Since there can be many attempts per call,
    * this method takes the last call attempt, and filters the metadata
    * key value pairs by those of interest.
    */
  private def convertMetadataToTestMap(metadataJson: String): Map[String, JsObject] = {
    val calls = metadataJson.parseJson.asJsObject.fields.get("calls").get.asJsObject
    val listOfkeys = Seq("runtimeAttributes", "preemptible", "executionStatus")

    val lastCallAttempt = calls.fields map {
      case (k1, allCallAttempts: JsArray) =>
        val lastCallAttempt = allCallAttempts.elements.last.asJsObject
        val filteredCallAttributes = lastCallAttempt.fields collect {
          case (k2, callAttributes) if listOfkeys.contains(k2) => k2 -> callAttributes
        }
        k1 -> JsObject(filteredCallAttributes)

      case(k, v) => throw new Exception(s"Wrong Js parameter for $k, Call attempts are expected to be in a JsArray")
    }

    lastCallAttempt
  }

  private def makeOutputMap(metadataJson: String): Map[String, JsValue] = {
    // Need DefaultJsonProtocol export to use convertTo[Map[String,String]]
    import DefaultJsonProtocol._
    val outputs = metadataJson.parseJson.asJsObject.fields.get("outputs").get.asJsObject.toString

    outputs.parseJson.convertTo[Map[String, JsValue]]
  }


  def verifyMetadataAndOutputs(workflow: Workflow, request: WorkflowRequest): Test[Workflow] = {
    new Test[Workflow] {
      val expectedMap: Option[Map[String, JsObject]] = request.metadata map { metadataString: String => convertMetadataToTestMap(metadataString) }
      val expectedOutputMap: Option[Map[String, JsValue]] = request.metadata map { metadataString: String => makeOutputMap(metadataString) }

      def verifyWorkflowMetadata(metadata: Map[String, JsObject]) = {
        expectedMap match {
          case Some(expected) if !expected.equals(metadata) => throw new Exception(s"Bad Metadata for Workflow ${request.name}")
          case _ =>
        }
      }

      def verifyWorkflowOutputs(outputs: Map[String, JsValue]) = {
        expectedOutputMap match {
          case Some(expected) if !expected.equals(outputs) => throw new Exception(s"Bad outputs for Workflow ${request.name}")
          case _ =>
        }
      }

      override def run: Try[Workflow] = {
        def ensureExpectedFile(metadata: String): Unit = {
          val expectedFile = request.base.getParent.resolve(s"${request.name}.metadata")
          if (!Files.exists(expectedFile)) {
            expectedFile.write(metadata)
          }
        }

        val response = MetadataRequest(Get(CentaurConfig.cromwellUrl + "/api/workflows/v1/" + workflow.id + "/metadata"))
        sendReceiveFutureCompletion(response map { allMetadata =>
          ensureExpectedFile(allMetadata)
          verifyWorkflowMetadata(convertMetadataToTestMap(allMetadata))
          verifyWorkflowOutputs(makeOutputMap(allMetadata))

          workflow
        })

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
