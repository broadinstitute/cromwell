package cromwell.engine.workflow.lifecycle.execution

import java.time.OffsetDateTime
import java.time.temporal.ChronoUnit
import java.util.concurrent.atomic.AtomicInteger

import akka.actor.ActorRef
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend._
import cromwell.core.ExecutionStatus._
import cromwell.core.io.AsyncIo
import cromwell.core.{JobKey, _}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData.DataStoreUpdate
import cromwell.engine.workflow.lifecycle.execution.keys._
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore.ValueKey
import cromwell.engine.workflow.lifecycle.execution.stores.{ExecutionStore, ValueStore}
import cromwell.engine.{EngineIoFunctions, EngineWorkflowDescriptor}
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

import scala.concurrent.ExecutionContext

object WorkflowExecutionDiff {
  def empty = WorkflowExecutionDiff(Map.empty)
}
/** Data differential between current execution data, and updates performed in a method that needs to be merged. */
final case class WorkflowExecutionDiff(executionStoreChanges: Map[JobKey, ExecutionStatus],
                                       jobKeyActorMappings: Map[ActorRef, JobKey] = Map.empty,
                                       valueStoreAdditions: Map[ValueKey, WomValue] = Map.empty,
                                       cumulativeOutputsChanges: Set[WomValue] = Set.empty,
                                       rootAndSubworkflowIds: Set[WorkflowId] = Set.empty) {
  def containsNewEntry: Boolean = {
    executionStoreChanges.exists(esc => esc._2 == NotStarted) || valueStoreAdditions.nonEmpty
  }
}

object WorkflowExecutionActorData {
  def apply(workflowDescriptor: EngineWorkflowDescriptor, ec: ExecutionContext, asyncIo: AsyncIo, totalJobsByRootWf: AtomicInteger, executionStore: ExecutionStore): WorkflowExecutionActorData = {
    WorkflowExecutionActorData(
      workflowDescriptor,
      executionStore,
      ValueStore.initialize(workflowDescriptor.knownValues),
      asyncIo,
      ec,
      totalJobsByRootWf = totalJobsByRootWf,
      rootAndSubworkflowIds = Set(workflowDescriptor.id)
    )
  }

  final case class DataStoreUpdate(runnableKeys: List[JobKey], statusChanges: Map[JobKey, ExecutionStatus], newData: WorkflowExecutionActorData)
}

case class WorkflowExecutionActorData(workflowDescriptor: EngineWorkflowDescriptor,
                                      executionStore: ExecutionStore,
                                      valueStore: ValueStore,
                                      asyncIo: AsyncIo,
                                      ec: ExecutionContext,
                                      jobKeyActorMappings: Map[ActorRef, JobKey] = Map.empty,
                                      jobFailures: Map[JobKey, Throwable] = Map.empty,
                                      downstreamExecutionMap: JobExecutionMap = Map.empty,
                                      totalJobsByRootWf: AtomicInteger,
                                      cumulativeOutputs: Set[WomValue] = Set.empty,
                                      rootAndSubworkflowIds: Set[WorkflowId]) extends StrictLogging {

  val expressionLanguageFunctions = new EngineIoFunctions(workflowDescriptor.pathBuilders, asyncIo, ec)

  def sealExecutionStore: WorkflowExecutionActorData = this.copy(
    executionStore = executionStore.seal
  )

  def callExecutionSuccess(jobKey: JobKey, outputs: CallOutputs, cumulativeOutputs: Set[WomValue], rootAndSubworkflowIds: Set[WorkflowId]): WorkflowExecutionActorData = {
    mergeExecutionDiff(WorkflowExecutionDiff(
      executionStoreChanges = Map(jobKey -> Done),
      valueStoreAdditions = toValuesMap(jobKey, outputs),
      cumulativeOutputsChanges = cumulativeOutputs ++ outputs.outputs.values,
      rootAndSubworkflowIds = rootAndSubworkflowIds
    ))
  }

  final def expressionEvaluationSuccess(expressionKey: ExpressionKey, values: Map[OutputPort, WomValue]): WorkflowExecutionActorData = {
    val valueStoreAdditions = values.map({
      case (outputPort, value) => ValueKey(outputPort, expressionKey.index) -> value
    })
    mergeExecutionDiff(WorkflowExecutionDiff(
      executionStoreChanges = Map(expressionKey -> Done),
      valueStoreAdditions =valueStoreAdditions
    ))
  }

  def executionFailure(failedJobKey: JobKey, reason: Throwable, jobExecutionMap: JobExecutionMap): WorkflowExecutionActorData = {
    mergeExecutionDiff(WorkflowExecutionDiff(
      executionStoreChanges = Map(failedJobKey -> ExecutionStatus.Failed))
    ).addExecutions(jobExecutionMap)
    .copy(
      jobFailures = jobFailures + (failedJobKey -> reason)
    )
  }

  def executionFailed(jobKey: JobKey): WorkflowExecutionActorData = mergeExecutionDiff(WorkflowExecutionDiff(Map(jobKey -> ExecutionStatus.Failed)))

  /** Converts call outputs to a ValueStore entries */
  private def toValuesMap(jobKey: JobKey, outputs: CallOutputs): Map[ValueKey, WomValue] = {
    outputs.outputs.map({
      case (outputPort, jobOutput) => ValueKey(outputPort, jobKey.index) -> jobOutput
    })
  }

  def addExecutions(jobExecutionMap: JobExecutionMap): WorkflowExecutionActorData = {
    this.copy(downstreamExecutionMap = downstreamExecutionMap ++ jobExecutionMap)
  }

  def removeJobKeyActor(actorRef: ActorRef): WorkflowExecutionActorData = {
    this.copy(
      jobKeyActorMappings = jobKeyActorMappings - actorRef
    )
  }

  def mergeExecutionDiff(diff: WorkflowExecutionDiff): WorkflowExecutionActorData = {
    val tsBefore = OffsetDateTime.now()
    val result = this.copy(
      executionStore = executionStore.updateKeys(diff.executionStoreChanges),
      valueStore = valueStore.add(diff.valueStoreAdditions),
      jobKeyActorMappings = jobKeyActorMappings ++ diff.jobKeyActorMappings,
      cumulativeOutputs = cumulativeOutputs ++ diff.cumulativeOutputsChanges,
      rootAndSubworkflowIds = rootAndSubworkflowIds ++ diff.rootAndSubworkflowIds
    )
    if (false) logger.warn(s"Merging $diff took ${ChronoUnit.MILLIS.between(OffsetDateTime.now(), tsBefore).abs} millis...")
    result
  }

  def mergeExecutionDiffs(diffs: Traversable[WorkflowExecutionDiff]): WorkflowExecutionActorData = {
    val tsBefore = OffsetDateTime.now()
    val result = preconsolidateDiffs(diffs) match {
      case Some(preconsolidated) =>
        logger.warn(s"Preconsolidation of ${diffs.size} diffs successful!")
        this.mergeExecutionDiff(preconsolidated)
      case None =>
        logger.warn(s"Unable to preconsolidation ${diffs.size} diffs")
        diffs.foldLeft(this)((newData, diff) => newData.mergeExecutionDiff(diff))
    }
    logger.warn(s"Merging ALL execution diffs took ${ChronoUnit.MILLIS.between(OffsetDateTime.now(), tsBefore).abs} millis...")
    result
  }

  def preconsolidateDiffs(diffs: Traversable[WorkflowExecutionDiff]): Option[WorkflowExecutionDiff] = {
    val result = diffs.foldLeft(WorkflowExecutionDiff.empty)((acc, diff) => {
      if (acc.executionStoreChanges.keySet.intersect(diff.executionStoreChanges.keySet).nonEmpty
        && acc.jobKeyActorMappings.keySet.intersect(diff.jobKeyActorMappings.keySet).nonEmpty
        && acc.valueStoreAdditions.keySet.intersect(diff.valueStoreAdditions.keySet).nonEmpty
        && acc.cumulativeOutputsChanges.intersect(diff.cumulativeOutputsChanges).nonEmpty
        && acc.rootAndSubworkflowIds.intersect(diff.rootAndSubworkflowIds).nonEmpty
      ) return None
      else {
        acc.copy(
          executionStoreChanges = acc.executionStoreChanges ++ diff.executionStoreChanges,
          jobKeyActorMappings = acc.jobKeyActorMappings ++ diff.jobKeyActorMappings,
          valueStoreAdditions = acc.valueStoreAdditions ++ diff.valueStoreAdditions,
          cumulativeOutputsChanges = acc.cumulativeOutputsChanges ++ diff.cumulativeOutputsChanges,
          rootAndSubworkflowIds = acc.rootAndSubworkflowIds ++ diff.rootAndSubworkflowIds
        )
      }
    })
    Option(result)
  }

  def jobExecutionMap: JobExecutionMap = {
    downstreamExecutionMap updated (workflowDescriptor.backendDescriptor, executionStore.startedJobs)
  }
  
  def executionStoreUpdate: DataStoreUpdate = {
    val update = executionStore.update
    DataStoreUpdate(update.runnableKeys, update.statusChanges, this.copy(executionStore = update.updatedStore))
  }

  def done: Boolean = executionStore.isDone
  def stalled: Boolean = executionStore.isStalled
}
