package wdl.transforms.draft2.wdlom2wom

import common.validation.ErrorOr.ErrorOr
import wdl.draft2.model.{WdlCallable, WdlTask, WdlWorkflow}
import wdl.draft2.model.{WdlTask, WdlWorkflow}
import wom.transforms.WomCallableMaker
import wom.callable.Callable
import wom.transforms.WomWorkflowDefinitionMaker.ops._
import wom.transforms.WomCommandTaskDefinitionMaker.ops._

object WdlDraft2WomCallableMaker extends WomCallableMaker[WdlCallable] {
  override def toWomCallable(callable: WdlCallable): ErrorOr[Callable] = callable match {
    case wf: WdlWorkflow => wf.toWomWorkflowDefinition(isASubworkflow = true)
    case t: WdlTask => t.toWomTaskDefinition
  }
}
