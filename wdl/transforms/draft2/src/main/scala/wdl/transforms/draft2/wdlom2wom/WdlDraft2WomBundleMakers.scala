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

object WdlDraft2WomBundleMakers {

  implicit val wdlDraft2NamespaceWomBundleMaker: WomBundleMaker[WdlNamespace] = new WomBundleMaker[WdlNamespace] {
    override def toWomBundle(from: WdlNamespace): Checked[WomBundle] = {
      val workflowsValidation: ErrorOr[List[WorkflowDefinition]] = from.workflows.toList.traverse(_.toWomWorkflowDefinition(isASubworkflow = false))
      val tasksValidation: ErrorOr[List[TaskDefinition]] = from.tasks.toList.traverse(_.toWomTaskDefinition)

      val errorOr = (workflowsValidation, tasksValidation) mapN { (workflows, tasks) =>
        val primary =
          if (workflows.size == 1) {
            workflows.headOption
          } else None
        WomBundle(primary, (tasks ++ workflows).map(c => c.name -> c).toMap, Map.empty, from.resolvedImportRecords) }
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
