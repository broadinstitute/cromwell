package cromwell.engine

import java.io.File

import akka.pattern.ask
import akka.util.Timeout
import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.WorkflowManagerActor.{WorkflowOutputs, _}
import cromwell.server.WorkflowManagerSystem
import cromwell.util.FileUtil
import spray.json._

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.concurrent.{Await, Promise}
import scala.language.postfixOps
import scala.util.Try

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
          if status.contains(WorkflowSucceeded)
          outputs <- (workflowManagerActor ? WorkflowOutputs(id)).mapTo[Option[Map[FullyQualifiedName, WdlValue]]]
        } yield outputs

        Await.result(futureOutputs, 5 seconds).get
      case _ => throw new RuntimeException("Expecting a JSON object")
    }
  }
}
