package centaur.api

import java.util.concurrent.Executors

import akka.actor.ActorSystem
import akka.stream.{ActorMaterializer, ActorMaterializerSettings}
import centaur.CentaurConfig
import centaur.test.metadata.WorkflowMetadata
import centaur.test.workflow.Workflow
import cromwell.api.CromwellClient
import cromwell.api.model.{CromwellBackends, SubmittedWorkflow, WorkflowStatus}

import scala.concurrent.{Await, ExecutionContext, Future}
import scala.concurrent.duration.FiniteDuration
import scala.util.Try

object CentaurCromwellClient {
  // Do not use scala.concurrent.ExecutionContext.Implicits.global as long as this is using Await.result
  // See https://github.com/akka/akka-http/issues/602
  // And https://github.com/viktorklang/blog/blob/master/Futures-in-Scala-2.12-part-7.md
  final implicit val blockingEc = ExecutionContext.fromExecutor(Executors.newCachedThreadPool())
  
  // Akka HTTP needs both the actor system and a materializer
  final implicit val system = ActorSystem("centaur-acting-like-a-system")
  final implicit val materializer: ActorMaterializer = ActorMaterializer(ActorMaterializerSettings(system))
  final val apiVersion = "v1"
  val cromwellClient = new CromwellClient(CentaurConfig.cromwellUrl, apiVersion)

  def submit(workflow: Workflow): Try[SubmittedWorkflow] = {
    sendReceiveFutureCompletion(cromwellClient.submit(workflow.toWorkflowSubmission(refreshToken = CentaurConfig.optionalToken)))
  }

  def status(workflow: SubmittedWorkflow): Try[WorkflowStatus] = {
    sendReceiveFutureCompletion(cromwellClient.status(workflow.id))
  }

  def metadata(workflow: SubmittedWorkflow): Try[WorkflowMetadata] = {
    sendReceiveFutureCompletion(cromwellClient.metadata(workflow.id)) map { m =>
      WorkflowMetadata.fromMetadataJson(m.value).toOption.get
    }
  }

  lazy val backends: Try[CromwellBackends] = sendReceiveFutureCompletion(cromwellClient.backends)

  /**
    * Ensure that the Future completes within the specified timeout. If it does not, or if the Future fails,
    * will return a Failure, otherwise a Success
    */
  def awaitFutureCompletion[T](x: Future[T], timeout: FiniteDuration) = Try(Await.result(x, timeout))
  def sendReceiveFutureCompletion[T](x: Future[T]) = awaitFutureCompletion(x, CentaurConfig.sendReceiveTimeout)
}
