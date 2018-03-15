package wdl.transforms.draft2.wdlom2wom

import cats.syntax.apply._
import cats.syntax.traverse._
import cats.instances.list._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import wdl.draft2.model.{WdlNamespace, WdlNamespaceWithWorkflow, WdlNamespaceWithoutWorkflow}
import wom.callable.{TaskDefinition, WorkflowDefinition}
import wom.executable.WomBundle
import wom.transforms.WomWorkflowDefinitionMaker.ops._
import wom.transforms.WomCommandTaskDefinitionMaker.ops._
import wom.transforms.WomBundleMaker

import scala.concurrent.Future

object WdlDraft2WomBundleMakers {

  implicit val wdlDraft2NamespaceWomBundleMaker: WomBundleMaker[WdlNamespace] = new WomBundleMaker[WdlNamespace] {
    override def toWomBundle(from: WdlNamespace): Checked[WomBundle] = {
      val workflowaValidation: ErrorOr[List[WorkflowDefinition]] = from.workflows.toList.traverse[ErrorOr, WorkflowDefinition](_.toWomWorkflowDefinition)
      val callsValidation: ErrorOr[List[TaskDefinition]] = from.tasks.toList.traverse[ErrorOr, TaskDefinition](_.toWomTaskDefinition)

      val errorOr = (workflowaValidation, callsValidation) mapN { (workflows, calls) => WomBundle((calls ++ workflows).toSet, Map.empty) }
      errorOr.toEither
    }
  }

  implicit val wdlDraft2NamespaceWithWorkflowWomBundleMaker: WomBundleMaker[WdlNamespaceWithWorkflow] = new WomBundleMaker[WdlNamespaceWithWorkflow] {
    override def toWomBundle(a: WdlNamespaceWithWorkflow): Checked[WomBundle] = wdlDraft2NamespaceWomBundleMaker.toWomBundle(a)
  }

  implicit val wdlDraft2NamespaceWithoutWorkflowWomBundleMaker: WomBundleMaker[WdlNamespaceWithoutWorkflow] = new WomBundleMaker[WdlNamespaceWithoutWorkflow] {
    override def toWomBundle(a: WdlNamespaceWithoutWorkflow): Checked[WomBundle] = wdlDraft2NamespaceWomBundleMaker.toWomBundle(a)
  }
}
