package cromwell.engine.workflow.lifecycle.execution

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.core.ExecutionStatus._
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.OutputStore.{OutputCallKey, OutputEntry}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{DeclarationKey, SubWorkflowKey}
import cromwell.engine.{EngineWorkflowDescriptor, WdlFunctions}
import wdl4s.values.WdlValue
import wdl4s.{GraphNode, Scope}

object WorkflowExecutionDiff {
  def empty = WorkflowExecutionDiff(Map.empty)
}
/** Data differential between current execution data, and updates performed in a method that needs to be merged. */
final case class WorkflowExecutionDiff(executionStoreChanges: Map[JobKey, ExecutionStatus],
                                       engineJobExecutionActorAdditions: Map[ActorRef, JobKey] = Map.empty) {
  def containsNewEntry = executionStoreChanges.exists(esc => esc._2 == NotStarted)
}

object WorkflowExecutionActorData {
  def empty(workflowDescriptor: EngineWorkflowDescriptor) = {
    new WorkflowExecutionActorData(
      workflowDescriptor,
      ExecutionStore.empty,
      Map.empty,
      Map.empty,
      Map.empty,
      Map.empty,
      OutputStore.empty
    )
  }
}

case class WorkflowExecutionActorData(workflowDescriptor: EngineWorkflowDescriptor,
                                      executionStore: ExecutionStore,
                                      backendJobExecutionActors: Map[JobKey, ActorRef],
                                      engineCallExecutionActors: Map[ActorRef, JobKey],
                                      subWorkflowExecutionActors: Map[SubWorkflowKey, ActorRef],
                                      downstreamExecutionMap: JobExecutionMap,
                                      outputStore: OutputStore) {

  val expressionLanguageFunctions = new WdlFunctions(workflowDescriptor.pathBuilders)

  def callExecutionSuccess(jobKey: JobKey, outputs: CallOutputs) = {
    val (newJobExecutionActors, newSubWorkflowExecutionActors) = jobKey match {
      case jobKey: BackendJobDescriptorKey => (backendJobExecutionActors - jobKey, subWorkflowExecutionActors)
      case swKey: SubWorkflowKey => (backendJobExecutionActors, subWorkflowExecutionActors - swKey)
      case _ => (backendJobExecutionActors, subWorkflowExecutionActors)
    }

    this.copy(
      executionStore = executionStore.add(Map(jobKey -> Done)),
      backendJobExecutionActors = newJobExecutionActors,
      subWorkflowExecutionActors = newSubWorkflowExecutionActors,
      outputStore = outputStore.add(updateSymbolStoreEntry(jobKey, outputs))
    )
  }

  def declarationEvaluationSuccess(declarationKey: DeclarationKey, value: WdlValue) = {
    val outputStoreKey = OutputCallKey(declarationKey.scope, declarationKey.index)
    val outputStoreValue = OutputEntry(declarationKey.scope.unqualifiedName, value.wdlType, Option(value))
    this.copy(
      executionStore = executionStore.add(Map(declarationKey -> Done)),
      outputStore = outputStore.add(Map(outputStoreKey -> List(outputStoreValue)))
    )
  }

  def executionFailed(jobKey: JobKey) = mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))

  /** Add the outputs for the specified `JobKey` to the symbol cache. */
  private def updateSymbolStoreEntry(jobKey: JobKey, outputs: CallOutputs) = {
    val newOutputEntries = outputs map {
      case (name, value) => OutputEntry(name, value.wdlValue.wdlType, Option(value.wdlValue))
    }

    Map(OutputCallKey(jobKey.scope, jobKey.index) -> newOutputEntries.toList)
  }

  /** Checks if the workflow is completed by scanning through the executionStore.
    * If complete, this will return Some(finalStatus).  Otherwise, returns None */
  def workflowCompletionStatus: Option[ExecutionStatus] = {
    // `List`ify the `prerequisiteScopes` to avoid expensive hashing of `Scope`s when assembling the result.

    def upstreamFailed(scope: Scope): List[GraphNode] = scope match {
      case node: GraphNode => node.upstreamAncestry.toList filter { s =>
        executionStore.store.exists({ case (key, status) => status == Failed && key.scope == s })
      }
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

  def removeEngineJobExecutionActor(actorRef: ActorRef) = {
    this.copy(engineCallExecutionActors = engineCallExecutionActors - actorRef)
  }

  def addCallExecutionActor(jobKey: JobKey, actor: Option[ActorRef]): WorkflowExecutionActorData = actor match {
      case Some(actorRef) =>
        jobKey match {
          case jobKey: BackendJobDescriptorKey => this.copy(backendJobExecutionActors = backendJobExecutionActors + (jobKey -> actorRef))
          case swKey: SubWorkflowKey => this.copy(subWorkflowExecutionActors = subWorkflowExecutionActors + (swKey -> actorRef))
          case _ => this
        }
      case None => this
  }

  def removeCallExecutionActor(jobKey: JobKey): WorkflowExecutionActorData = {
    jobKey match {
      case jobKey: BackendJobDescriptorKey => this.copy(backendJobExecutionActors = backendJobExecutionActors - jobKey)
      case swKey: SubWorkflowKey => this.copy(subWorkflowExecutionActors = subWorkflowExecutionActors - swKey)
      case _ => this
    }
  }

  def addExecutions(jobExecutionMap: JobExecutionMap): WorkflowExecutionActorData = {
    this.copy(downstreamExecutionMap = downstreamExecutionMap ++ jobExecutionMap)
  }

  def mergeExecutionDiff(diff: WorkflowExecutionDiff): WorkflowExecutionActorData = {
    this.copy(
      executionStore = executionStore.add(diff.executionStoreChanges),
      engineCallExecutionActors = engineCallExecutionActors ++ diff.engineJobExecutionActorAdditions)
  }

  def mergeExecutionDiffs(diffs: Traversable[WorkflowExecutionDiff]): WorkflowExecutionActorData = {
    diffs.foldLeft(this)((newData, diff) => newData.mergeExecutionDiff(diff))
  }
  
  def jobExecutionMap: JobExecutionMap = {
    val keys = executionStore.store.collect({case (k: BackendJobDescriptorKey, status) if status != ExecutionStatus.NotStarted => k }).toList
    downstreamExecutionMap updated (workflowDescriptor.backendDescriptor, keys)
  }
  
  def hasRunningActors = backendJobExecutionActors.nonEmpty || subWorkflowExecutionActors.nonEmpty
}
