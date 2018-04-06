package wdl.transforms.draft2.wdlom2wom

// TODO: 2.11: cats.syntax.either._ AND typeclass longhand form, oh my!
import cats.syntax.either._
import common.Checked
import wdl.draft2.model.{WdlNamespace, WdlNamespaceWithWorkflow, WdlNamespaceWithoutWorkflow}
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomBundleMakers._
import wdl.shared.transforms.wdlom2wom.WdlSharedInputParsing
import wom.core.WorkflowJson
import wom.executable.Executable
import wom.expression.IoFunctionSet
import wom.transforms.WomExecutableMaker
import wom.transforms.WomBundleMaker.ops._

object WdlDraft2WomExecutableMakers {
  implicit val namespaceWomExecutableMaker: WomExecutableMaker[WdlNamespace] = new WomExecutableMaker[WdlNamespace] {
    override def toWomExecutable(a: WdlNamespace, inputs: Option[WorkflowJson], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] = {
      a.toWomBundle flatMap { bundle =>
        WdlSharedInputParsing.buildWomExecutable(bundle, inputs, ioFunctions, strictValidation)
      }
    }
  }

  implicit val namespaceWithWorkflowWomExecutableMaker: WomExecutableMaker[WdlNamespaceWithWorkflow] = new WomExecutableMaker[WdlNamespaceWithWorkflow] {
    override def toWomExecutable(a: WdlNamespaceWithWorkflow, inputs: Option[WorkflowJson], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] =
      namespaceWomExecutableMaker.toWomExecutable(a, inputs, ioFunctions, strictValidation)
  }

  implicit val namespaceWithoutWorkflowWomExecutableMaker: WomExecutableMaker[WdlNamespaceWithoutWorkflow] = new WomExecutableMaker[WdlNamespaceWithoutWorkflow] {
    override def toWomExecutable(a: WdlNamespaceWithoutWorkflow, inputs: Option[WorkflowJson], ioFunctions: IoFunctionSet, strictValidation: Boolean): Checked[Executable] =
      namespaceWomExecutableMaker.toWomExecutable(a, inputs, ioFunctions, strictValidation)
  }
}
