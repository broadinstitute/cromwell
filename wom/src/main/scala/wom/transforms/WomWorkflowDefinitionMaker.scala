package wom.transforms

import common.validation.ErrorOr.ErrorOr
import wom.callable.WorkflowDefinition

trait WomWorkflowDefinitionMaker[A] {
  def toWomWorkflowDefinition(a: A): ErrorOr[WorkflowDefinition]
}

object WomWorkflowDefinitionMaker {
  // This apply lets us grab an appropriate WomXMaker[A] out of implicit scope like "val maker = WomXMaker[A]"
  // eg used in the implicit class below.
  def apply[A](implicit maker: WomWorkflowDefinitionMaker[A]): WomWorkflowDefinitionMaker[A] = maker

  // The restriction [A: WomXMaker] is scala syntax magic for "if there exists in scope a WomXMaker for A"
  implicit class CanMakeWorkflowDefinition[A: WomWorkflowDefinitionMaker](val a: A) {
    def toWomWorkflowDefinition = WomWorkflowDefinitionMaker[A].toWomWorkflowDefinition(a)
  }
}
