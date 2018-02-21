package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.collections.EnhancedCollections._
import wdl.draft3.parser.WdlParser.Ast
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{InputsSectionElement, OutputsSectionElement, WorkflowBodyElement, WorkflowDefinitionElement}


object AstToWorkflowDefinitionElement {

  def convert(a: Ast): ErrorOr[WorkflowDefinitionElement] = {

    val nameElementValidation: ErrorOr[String] = astNodeToString(a.getAttribute("name")).toValidated

    val bodyElementsValidation: ErrorOr[Vector[WorkflowBodyElement]] = a.getAttributeAsVector[WorkflowBodyElement]("body")(astNodeToWorkflowBodyElement).toValidated

    (nameElementValidation, bodyElementsValidation) flatMapN combineElements
  }

  private def combineElements(name: String, bodyElements: Vector[WorkflowBodyElement]) = {

      val inputsSectionValidation: ErrorOr[Option[InputsSectionElement]] = validateSize(bodyElements.filterByType[InputsSectionElement], "inputs", 1)
      val outputsSectionValidation: ErrorOr[Option[OutputsSectionElement]] = validateSize(bodyElements.filterByType[OutputsSectionElement], "outputs", 1)

      (inputsSectionValidation, outputsSectionValidation) mapN { (validInputs, validOutputs) =>
        WorkflowDefinitionElement(name, validInputs, validOutputs)
      }
  }

  private def validateSize[A](elements: Vector[A], sectionName: String, numExpected: Int): ErrorOr[Option[A]] = {
    val sectionValidation: ErrorOr[Option[A]] = if (elements.size > numExpected) {
      s"Workflow cannot have more than $numExpected $sectionName sections, found ${elements.size}.".invalidNel
    } else {
      elements.headOption.validNel
    }
    sectionValidation
  }
}
