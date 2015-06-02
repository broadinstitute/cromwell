package cromwell.engine

import java.io.File

import cromwell.binding.FullyQualifiedName
import cromwell.binding.types._
import cromwell.binding.values.WdlValue
import cromwell.engine.WorkflowManagerActor._
import cromwell.server.WorkflowManagerSystem
import cromwell.util.{ActorUtil, FileUtil}
import spray.json._
import scala.concurrent.Future
import scala.util.{Try, Failure}

class SingleWorkflowRunner extends WorkflowManagerSystem {
  def run(wdlFile: File, inputs: File): Try[Map[FullyQualifiedName, WdlValue]] = {
    Try(FileUtil.slurp(inputs).parseJson) map {
      case JsObject(rawInputs) =>
        val workflowId = ActorUtil.messageAndWait(SubmitWorkflow(FileUtil.slurp(wdlFile), rawInputs), _.mapTo[WorkflowId])(workflowManagerActor)
        println(s"Workflow ID: $workflowId")

        /* TODO: horrible, terrible hack to wait until a workflow is succeeded */
        workflowManagerActor ! NotifyCompletion(this)
        this.synchronized {
          try {
            this.wait()
          } catch {
            case _: InterruptedException => println("interrupted")
          }
        }

        val state = ActorUtil.messageWaitAndGet(WorkflowStatus(workflowId), _.mapTo[Option[WorkflowState]])(workflowManagerActor)
        println(s"Workflow state: $state")
        state match {
          case WorkflowSucceeded =>
            val outputs = ActorUtil.messageWaitAndGet(WorkflowOutputs(workflowId), _.mapTo[Option[Map[FullyQualifiedName, WdlValue]]])(workflowManagerActor)
            workflowManagerActor ! Shutdown()
            outputs
          case _ =>
            throw new RuntimeException("Workflow Failed")
        }
      case _ => throw new RuntimeException("Expecting a JSON object")
    }
  }
}
