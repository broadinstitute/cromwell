package cromwell.engine.workflow.lifecycle.execution

import common.validation.ErrorOr.ErrorOr
import cromwell.core.{ExecutionStatus, JobKey}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore.ValueKey
import wom.values.WomOptionalValue

package object keys {

  def processBypassedNode(jobKey: JobKey, data: WorkflowExecutionActorData): ErrorOr[WorkflowExecutionDiff] = {
    import cats.syntax.validated._
    WorkflowExecutionDiff(
      executionStoreChanges = Map(jobKey -> ExecutionStatus.Bypassed),
      valueStoreAdditions = bypassedScopeResults(jobKey)
    ).validNel
  }

  private def bypassedScopeResults(jobKey: JobKey): Map[ValueKey, WomOptionalValue] =
    jobKey.node.outputPorts.map { outputPort =>
      ValueKey(outputPort, jobKey.index) -> WomOptionalValue.none(outputPort.womType)
    }.toMap
}
