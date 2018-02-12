package wdl.draft3.transforms.wdlom2wom

import cats.data.Validated.Valid
import org.scalatest.{FlatSpec, Matchers}
import wdl.model.draft3.elements.WorkflowDefinitionElement

class WorkflowDefinitionElementToWomWorkflowDefinitionSpec extends FlatSpec with Matchers {

  val cases = List(
    ("empty workflow definition", WorkflowDefinitionElement("empty_workflow"), None)
  )


  cases foreach { case (testName, workflowDefinitionElement, inputs) =>
    it should s"convert the '$testName' into WOM" in {

      val convertor = WorkflowDefinitionElementToWomWorkflowDefinition(inputs)
      convertor.convert(workflowDefinitionElement) match {
        case Valid(_) => // Great!
        case errors => fail(s"Failed to produce WOM for '$testName': ${errors.toList.mkString(System.lineSeparator, System.lineSeparator, System.lineSeparator)}")
      }
    }
  }
}
