package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.callable.WorkflowDefinition
import simulacrum._

@typeclass
trait WomWorkflowDefinitionMaker[A] {
  def toWomWorkflowDefinition(a: A, isASubworkflow: Boolean): ErrorOr[WorkflowDefinition]
}
