package cromwell.engine.workflow.lifecycle.execution

import _root_.wdl.values.WdlValue
import akka.actor.ActorRef
import cromwell.backend._
import cromwell.core.ExecutionStatus._
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.ValueStore.ValueKey
import cromwell.engine.workflow.lifecycle.execution.keys._
import cromwell.engine.{EngineWorkflowDescriptor, WdlFunctions}

import scala.language.postfixOps

object WorkflowExecutionDiff {
  def empty = WorkflowExecutionDiff(Map.empty)
}
/** Data differential between current execution data, and updates performed in a method that needs to be merged. */
final case class WorkflowExecutionDiff(executionStoreChanges: Map[JobKey, ExecutionStatus],
                                       engineJobExecutionActorAdditions: Map[ActorRef, JobKey] = Map.empty,
                                       valueStoreAdditions: Map[ValueKey, WdlValue] = Map.empty) {
  def containsNewEntry: Boolean = {
    executionStoreChanges.exists(esc => esc._2 == NotStarted) || valueStoreAdditions.nonEmpty
  }
}

object WorkflowExecutionActorData {
  def empty(workflowDescriptor: EngineWorkflowDescriptor): WorkflowExecutionActorData = {
    new WorkflowExecutionActorData(
      workflowDescriptor,
      ExecutionStore.empty,
      Map.empty,
      Map.empty,
      Map.empty,
      Map.empty,
      ValueStore.empty
    )
  }
}

case class WorkflowExecutionActorData(workflowDescriptor: EngineWorkflowDescriptor,
                                      executionStore: ExecutionStore,
                                      backendJobExecutionActors: Map[JobKey, ActorRef],
                                      engineCallExecutionActors: Map[ActorRef, JobKey],
                                      subWorkflowExecutionActors: Map[SubWorkflowKey, ActorRef],
                                      downstreamExecutionMap: JobExecutionMap,
                                      valueStore: ValueStore) {

  val expressionLanguageFunctions = new WdlFunctions(workflowDescriptor.pathBuilders)

  def isInBypassedScope(jobKey: JobKey): Boolean = {
    executionStore.isInBypassedConditional(jobKey)
  }

  def callExecutionSuccess(jobKey: JobKey, outputs: CallOutputs): WorkflowExecutionActorData = {
    val (newJobExecutionActors, newSubWorkflowExecutionActors) = jobKey match {
      case jobKey: BackendJobDescriptorKey => (backendJobExecutionActors - jobKey, subWorkflowExecutionActors)
      case swKey: SubWorkflowKey => (backendJobExecutionActors, subWorkflowExecutionActors - swKey)
      case _ => (backendJobExecutionActors, subWorkflowExecutionActors)
    }

    this.copy(
      executionStore = executionStore.add(Map(jobKey -> Done)),
      backendJobExecutionActors = newJobExecutionActors,
      subWorkflowExecutionActors = newSubWorkflowExecutionActors,
      valueStore = valueStore.add(updateSymbolStoreEntry(jobKey, outputs))
    )
  }

  final def expressionEvaluationSuccess(expressionKey: ExpressionKey, value: WdlValue): WorkflowExecutionActorData = {
    val valueStoreKey = ValueKey(expressionKey.singleOutputPort, expressionKey.index)
    this.copy(
      executionStore = executionStore.add(Map(expressionKey -> Done)),
      valueStore = valueStore.add(Map(valueStoreKey -> value))
    )
  }

  def executionFailed(jobKey: JobKey): WorkflowExecutionActorData = mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))

  /** Add the outputs for the specified `JobKey` to the symbol cache. */
  private def updateSymbolStoreEntry(jobKey: JobKey, outputs: CallOutputs): Map[ValueKey, WdlValue] = {
    jobKey.node.outputPorts flatMap { outputPort =>
      outputs.collectFirst { 
        case (name, JobOutput(value)) if name == outputPort.name => value
      } map { ValueKey(outputPort, jobKey.index) -> _ }
    } toMap
  }

  /** Checks if the workflow is completed.
    * If complete, this will return Some(finalStatus).  Otherwise, returns None */
  def workflowCompletionStatus: Option[ExecutionStatus] = {
    (executionStore.hasActiveJob, executionStore.hasFailedJob) match {
      case (false, true) => Option(Failed)
      case (false, false) => Option(Done)
      case _ => None
    }
  }

  def removeEngineJobExecutionActor(actorRef: ActorRef): WorkflowExecutionActorData = {
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
      valueStore = valueStore.add(diff.valueStoreAdditions),
      engineCallExecutionActors = engineCallExecutionActors ++ diff.engineJobExecutionActorAdditions
    )
  }

  def mergeExecutionDiffs(diffs: Traversable[WorkflowExecutionDiff]): WorkflowExecutionActorData = {
    diffs.foldLeft(this)((newData, diff) => newData.mergeExecutionDiff(diff))
  }
  
  def resetCheckRunnable: WorkflowExecutionActorData = this.copy(executionStore = executionStore.copy(hasNewRunnables = false))
  
  def hasNewRunnables: Boolean = executionStore.hasNewRunnables
  
  def jobExecutionMap: JobExecutionMap = {
    downstreamExecutionMap updated (workflowDescriptor.backendDescriptor, executionStore.startedJobs)
  }
  
  def hasRunningActors: Boolean = backendJobExecutionActors.nonEmpty || subWorkflowExecutionActors.nonEmpty
}
