package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.callable.WorkflowDefinition
import simulacrum._
import scala.language.implicitConversions

@typeclass
trait WomWorkflowDefinitionMaker[A] {
  @op("toWomWorkflowDefinition")
  def toWomWorkflowDefinition(a: A): ErrorOr[WorkflowDefinition]
}
