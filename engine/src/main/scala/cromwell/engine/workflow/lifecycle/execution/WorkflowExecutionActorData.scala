package cromwell.engine.workflow.lifecycle.execution

import akka.actor.ActorRef
import cromwell.backend.JobKey
import cromwell.core._
import cromwell.engine.ExecutionStatus._
import cromwell.engine.workflow.lifecycle.execution.OutputStore.{OutputCallKey, OutputEntry}
import cromwell.engine.{EngineWorkflowDescriptor, ExecutionStatus, WdlFunctions}
import cromwell.webservice.WdlValueJsonFormatter

import scala.language.postfixOps

object WorkflowExecutionDiff {
  def empty = WorkflowExecutionDiff(Map.empty)
}
/** Data differential between current execution data, and updates performed in a method that needs to be merged. */
final case class WorkflowExecutionDiff(executionStore: Map[JobKey, ExecutionStatus]) {
  def containsNewEntry = executionStore.exists(_._2 == NotStarted)
}

case class WorkflowExecutionActorData(workflowDescriptor: EngineWorkflowDescriptor,
                                      executionStore: ExecutionStore,
                                      backendJobExecutionActors: Map[JobKey, ActorRef],
                                      outputStore: OutputStore) extends WdlLookup {

  override val expressionLanguageFunctions = new WdlFunctions(workflowDescriptor.backendDescriptor.workflowOptions)

  def jobExecutionSuccess(jobKey: JobKey, outputs: JobOutputs) = this.copy(
    executionStore = executionStore.add(Map(jobKey -> Done)),
    backendJobExecutionActors = backendJobExecutionActors - jobKey,
    outputStore = outputStore.add(updateSymbolStoreEntry(jobKey, outputs))
  )

  /** Add the outputs for the specified `JobKey` to the symbol cache. */
  private def updateSymbolStoreEntry(jobKey: JobKey, outputs: JobOutputs) = {
    val newOutputEntries = outputs map {
      case (name, value) => OutputEntry(name, value.wdlValue.wdlType, Option(value.wdlValue))
    }

    Map(OutputCallKey(jobKey.scope, jobKey.index) -> newOutputEntries)
  }

  /** Checks if the workflow is completed by scanning through the executionStore */
  def isWorkflowComplete: Boolean = {
    def isDone(executionStatus: ExecutionStatus): Boolean = executionStatus == Done || executionStatus == Preempted
    executionStore.store.values.forall(isDone)
  }

  def hasFailedJob: Boolean = {
    executionStore.store.values.exists(_ == ExecutionStatus.Failed)
  }

  def addBackendJobExecutionActor(jobKey: JobKey, actor: ActorRef): WorkflowExecutionActorData = {
    this.copy(backendJobExecutionActors = backendJobExecutionActors + (jobKey -> actor))
  }

  def removeBackendJobExecutionActor(jobKey: JobKey): WorkflowExecutionActorData = {
    this.copy(backendJobExecutionActors = backendJobExecutionActors - jobKey)
  }

  def outputsJson(): String = {
    // Printing the final outputs, temporarily here until SingleWorkflowManagerActor is made in-sync with the shadow mode
    import WdlValueJsonFormatter._
    import spray.json._
    val workflowOutputs = outputStore.store collect {
      case (key, outputs) if key.index.isEmpty => outputs map { output =>
        s"${key.call.fullyQualifiedName}.${output.name}" -> (output.wdlValue map { _.valueString } getOrElse "N/A")
      }
    }

    "Workflow complete. Final Outputs: \n" + workflowOutputs.flatten.toMap.toJson.prettyPrint
  }

  def mergeExecutionDiff(diff: WorkflowExecutionDiff): WorkflowExecutionActorData = {
    this.copy(executionStore = executionStore.add(diff.executionStore))
  }

  def mergeExecutionDiffs(diffs: Traversable[WorkflowExecutionDiff]): WorkflowExecutionActorData = {
    diffs.foldLeft(this)((newData, diff) => newData.mergeExecutionDiff(diff))
  }

}
