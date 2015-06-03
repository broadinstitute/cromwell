package cromwell.engine

import java.io.File

import akka.actor.{ActorRef, Actor}
import akka.event.{LoggingReceive, Logging}
import akka.pattern.ask
import akka.util.Timeout
import cromwell.binding._
import cromwell.binding.types._
import cromwell.binding.values.WdlValue
import cromwell.engine.StoreActor.UpdateStatus
import cromwell.engine.WorkflowManagerActor.WorkflowOutputs
import cromwell.engine.WorkflowManagerActor._
import cromwell.engine.backend.Backend
import cromwell.server.WorkflowManagerSystem
import cromwell.util.{ActorUtil, FileUtil}
import spray.json._
import scala.concurrent.{Await, Promise, ExecutionContext, Future}
import scala.util.{Success, Try, Failure}
import scala.concurrent.duration._
import ExecutionContext.Implicits.global
import scala.language.postfixOps

class SingleWorkflowRunner extends WorkflowManagerSystem {
  implicit val timeout = Timeout(5 seconds)
  def run(wdlFile: File, inputs: File): Try[Map[FullyQualifiedName, WdlValue]] = {
    Try(FileUtil.slurp(inputs).parseJson) map {
      case JsObject(rawInputs) =>
        val futureWorkflowId = (workflowManagerActor ? SubmitWorkflow(FileUtil.slurp(wdlFile), rawInputs)).mapTo[WorkflowId]
        val id = Await.result(futureWorkflowId, 5 seconds)
        val promise = Promise[Unit]()
        println(s"Workflow ID: $id")
        workflowManagerActor ? NotifyCompletion(promise)
        Await.result(promise.future, 5 seconds)

        val futureOutputs = for {
          status <- (workflowManagerActor ? WorkflowStatus(id)).mapTo[Option[WorkflowState]]
          if status == Some(WorkflowSucceeded)
          outputs <- (workflowManagerActor ? WorkflowOutputs(id)).mapTo[Option[Map[FullyQualifiedName, WdlValue]]]
        } yield outputs

        Await.result(futureOutputs, 5 seconds).get
      case _ => throw new RuntimeException("Expecting a JSON object")
    }
  }
}
