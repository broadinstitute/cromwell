package wdl.draft3.transforms.wdlom2wom

import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements.WorkflowDefinitionElement

class WorkflowDefinitionElementToWomWorkflowDefinitionSpec extends FlatSpec with Matchers {

  val cases = List(
    ("empty workflow definition", WorkflowDefinitionElement("empty_workflow"))
  )


  cases foreach { case (testName, workflowDefinitionElement) =>
    it should s"convert the '$testName' into WOM" in {

      checkedWorkflowDefinitionElementToWomWorkflowDefinition.run(workflowDefinitionElement) match {
        case Right(_) => // Great!
        case Left(errors) => fail(s"Failed to produce WOM for '$testName': ${errors.toList.mkString(System.lineSeparator, System.lineSeparator, System.lineSeparator)}")
      }
    }
  }
}
