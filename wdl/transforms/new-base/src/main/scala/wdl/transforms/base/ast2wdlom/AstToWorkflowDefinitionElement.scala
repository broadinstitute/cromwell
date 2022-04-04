package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import common.collections.EnhancedCollections._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import wdl.model.draft3.elements._
import wom.SourceFileLocation
import wdl.model.draft3.elements.ExpressionElement.StdoutElement

object AstToWorkflowDefinitionElement {

  def astToWorkflowDefinitionElement(implicit astNodeToWorkflowBodyElement: CheckedAtoB[GenericAstNode, WorkflowBodyElement]
                                    ): CheckedAtoB[GenericAst, WorkflowDefinitionElement] = CheckedAtoB.fromErrorOr { a: GenericAst =>
    val nameElementValidation: ErrorOr[String] = astNodeToString(a.getAttribute("name")).toValidated

    val sourceLocation : Option[SourceFileLocation] = a.getSourceLine.map(SourceFileLocation(_))
    val bodyElementsValidation: ErrorOr[Vector[WorkflowBodyElement]] = a.getAttributeAsVector[WorkflowBodyElement]("body")(astNodeToWorkflowBodyElement).toValidated
    (nameElementValidation, bodyElementsValidation) flatMapN { (name, bodyElements) =>
      combineElements(name, sourceLocation, bodyElements)
    }
  }

  private def combineElements(name: String,
                              sourceLocation: Option[SourceFileLocation],
                              bodyElements: Vector[WorkflowBodyElement]) = {

    val inputsSectionValidation: ErrorOr[Option[InputsSectionElement]] = validateSize(bodyElements.filterByType[InputsSectionElement], "inputs", 1)
    val outputsSectionValidation: ErrorOr[Option[OutputsSectionElement]] = validateSize(bodyElements.filterByType[OutputsSectionElement], "outputs", 1)

    val inputsSectionStdoutValidation: ErrorOr[ErrorOr[Option[String]]] = inputsSectionValidation.map(checkIfStdoutExist)

    val graphSections: Vector[WorkflowGraphElement] = bodyElements.filterByType[WorkflowGraphElement]

    val metaSectionValidation: ErrorOr[Option[MetaSectionElement]] = validateSize(bodyElements.filterByType[MetaSectionElement], "meta", 1)
    val parameterMetaSectionValidation: ErrorOr[Option[ParameterMetaSectionElement]] = validateSize(bodyElements.filterByType[ParameterMetaSectionElement], "parameterMeta", 1)

    (inputsSectionValidation, outputsSectionValidation, metaSectionValidation, parameterMetaSectionValidation, inputsSectionStdoutValidation) mapN {
      (validInputs, validOutputs, meta, parameterMeta, _) =>
      WorkflowDefinitionElement(name, validInputs, graphSections.toSet, validOutputs, meta, parameterMeta, sourceLocation)
    }
  }

  private def checkIfStdoutExist(inputSection: Option[InputsSectionElement]): ErrorOr[Option[String]] = {
    inputSection match {
      case Some(section) =>
        if (section.inputDeclarations.flatMap(_.expression).exists(_.isInstanceOf[StdoutElement.type])) {
          s"Workflow cannot have stdout expression in input section at workflow-level.".invalidNel
        } else None.validNel
      case None => None.validNel
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

  /*for(validate<-AstToWorkflowDefinitionElement){
    validate match {
      case check input{}
      case check p

    }
    */

}

//def validateStdout(a:)
