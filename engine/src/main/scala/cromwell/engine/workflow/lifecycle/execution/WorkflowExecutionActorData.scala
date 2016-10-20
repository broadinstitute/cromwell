package cromwell.engine.workflow.lifecycle.execution

import akka.actor.ActorRef
import cromwell.core.ExecutionStatus._
import cromwell.core.OutputStore.{OutputCallKey, OutputEntry}
import cromwell.core._
import cromwell.engine.{EngineWorkflowDescriptor, WdlFunctions}
import cromwell.util.JsonFormatting.WdlValueJsonFormatter
import wdl4s.Scope


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

  override val expressionLanguageFunctions = new WdlFunctions(workflowDescriptor.engineFilesystems)

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

  /** Checks if the workflow is completed by scanning through the executionStore.
    * If complete, this will return Some(finalStatus).  Otherwise, returns None */
  def workflowCompletionStatus: Option[ExecutionStatus] = {
    // `List`ify the `prerequisiteScopes` to avoid expensive hashing of `Scope`s when assembling the result.
    def upstream(scope: Scope): List[Scope] = scope.prerequisiteScopes.toList ++ scope.prerequisiteScopes.toList.flatMap(upstream)
    def upstreamFailed(scope: Scope) = upstream(scope) filter { s =>
      executionStore.store.contains({ case (key: JobKey, status: ExecutionStatus) => status == Failed && key.scope == s })
    }
    // activeJobs is the subset of the executionStore that are either running or will run in the future.
    val activeJobs = executionStore.store.toList filter {
      case (jobKey, jobStatus) => (jobStatus == NotStarted && upstreamFailed(jobKey.scope).isEmpty) || jobStatus == QueuedInCromwell || jobStatus == Starting || jobStatus == Running
    }

    activeJobs match {
      case jobs if jobs.isEmpty && hasFailedJob => Option(Failed)
      case jobs if jobs.isEmpty && !hasFailedJob => Option(Done)
      case _ => None
    }
  }

  def hasFailedJob: Boolean = {
    executionStore.store.values.exists(_ == ExecutionStatus.Failed)
  }

  def addBackendJobExecutionActor(jobKey: JobKey, actor: Option[ActorRef]): WorkflowExecutionActorData = actor match {
      case Some(actorRef) => this.copy(backendJobExecutionActors = backendJobExecutionActors + (jobKey -> actorRef))
      case None => this
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
