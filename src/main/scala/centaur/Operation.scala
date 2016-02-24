package centaur

import java.util.UUID

import akka.actor.ActorSystem
import spray.client.pipelining._
import spray.http.{HttpRequest, FormData}
import spray.httpx.SprayJsonSupport._

import scala.annotation.tailrec
import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{Future, Await}
import scala.util.{Failure, Success, Try}
import scala.concurrent.ExecutionContext.Implicits.global

import CromwellStatusJsonSupport._
import Operation._

/**
  * Basic building blocks which can be chained together to form a TestFormula. The idea here is to provide a mini-DSL
  * which will allow someone to chain these operations together as part of a new test plan. As an example:
  * "submit a workflow, wait until it is complete, verify its outputs". These components can be chained together by
  * a Formula via a for comprehension
  */
sealed trait Operation {
  /**
    * FIXME: Right now this scheme is piggybacking on Try and so I need this perform here. With a different monad
    * structure we could do away with the perform but i wanted to get the pieces in place ASAP
    */
  def perform: Test
}

/**
  * Submits a workflow to Cromwell and ensures that the ID was received within the requested timeout
  */
final case class SubmitWorkflow(request: WorkflowRequest) extends Operation {
  override def perform: Test = {
    val formData = FormData(List("wdlSource" -> request.wdl, "workflowInputs" -> request.inputs, "workflowOptions" -> request.options))
    val response = Pipeline(Post(CentaurConfig.cromwellUrl + "/api/workflows/v1", formData))
    sendReceiveFutureCompletion(response map { _.id } map UUID.fromString map { Workflow(_, CentaurConfig.cromwellUrl) })
  }
}

final case class PollUntilStatus(workflow: Workflow, expectedStatus: WorkflowStatus) extends Operation {
  def perform: Test = workflowLengthFutureCompletion(Future { doPerform() })

  @tailrec
  private def doPerform(): Workflow = {
    val response = Pipeline(Get(CentaurConfig.cromwellUrl + "/api/workflows/v1/" + workflow.id + "/status"))
    val status = sendReceiveFutureCompletion(response map { r => WorkflowStatus(r.status) })
    status match {
      case Success(s) if s == expectedStatus => workflow
      case Failure(f) => throw f
      case _ =>
        Thread.sleep(10000) // This could be a lot smarter including cromwell style backoff
        doPerform()
    }
  }
}

final case class PollUntilTerminal(workflow: Workflow, terminationStatus: TerminalStatus) extends Operation {
  override def perform: Test =  PollUntilStatus(workflow, terminationStatus).perform
}

/*
 TODO: Example operations
  case class ValidateOutputs(workflow: Workflow, outputs: SomePlaceholder) extends Operation
  case class ValidateMetadata(workflow: Workflow, metadata: SomePlaceholder) extends Operation
  case class ValidateFileExists(workflow: Workflow, filePath: SomePlaceholderWhichcouldBeLocal/GCS/Etc) extends Operation
  case class Delay(timePlaceholder: Int)
*/

object Operation {
  /*
   * Monad wrapping a workflow under test. Currently piggybacking on Try although something which allows chaining
   * w/o needing to call an extra method would be ideal
   *
   * FIXME/TODO - Miguel pointed out that this would be better using a Future where operations which are actually
   * using Futures would use the 'after' pattern internally instead of the awaits below. That's true and perhaps
   * how things will shake out but we decided to hold off until the final structure of the operation chaining is
   * known
   */
  type Test = Try[Workflow]

  // Spray needs an implicit ActorSystem
  implicit val system = ActorSystem("centaur")

  // FIXME: Pretty sure this will be insufficient once we move past submit & polling, but hey, continuous improvement!
  val Pipeline: HttpRequest => Future[CromwellStatus] = sendReceive ~> unmarshal[CromwellStatus]

  /**
    * Ensure that the Future completes within the specified timeout. If it does not, or if the Future fails,
    * will return a Failure, otherwise a Success
    */
  def awaitFutureCompletion[T](x: Future[T], timeout: FiniteDuration) = Try(Await.result(x, timeout))
  def sendReceiveFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.sendReceiveTimeout)
  def workflowLengthFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.maxWorkflowLength)
}

