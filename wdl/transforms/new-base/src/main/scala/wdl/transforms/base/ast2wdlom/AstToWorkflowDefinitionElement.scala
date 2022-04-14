package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import common.collections.EnhancedCollections._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import wdl.model.draft3.elements._
import wom.SourceFileLocation
import wdl.model.draft3.elements.ExpressionElement._

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

    val inputsSectionValidation: ErrorOr[Option[InputsSectionElement]] = for {
      inputValidateElement <- validateSize(bodyElements.filterByType[InputsSectionElement], "inputs", 1): ErrorOr[Option[InputsSectionElement]]
      _ <- checkIfStdoutInputExist(inputValidateElement)
      _ <- checkIfStderrInputExist(inputValidateElement)
    } yield inputValidateElement

    val intermediateValueDeclarationStdoutCheck: ErrorOr[Option[String]] = checkStdoutIntermediates(bodyElements.filterByType[IntermediateValueDeclarationElement])
    val intermediateValueDeclarationStderrCheck: ErrorOr[Option[String]] = checkStderrIntermediates(bodyElements.filterByType[IntermediateValueDeclarationElement])

    val outputsSectionValidation: ErrorOr[Option[OutputsSectionElement]] = for {
      outputValidateElement <- validateSize(bodyElements.filterByType[OutputsSectionElement], "outputs", 1): ErrorOr[Option[OutputsSectionElement]]
      _ <- checkIfStdoutOutputExists(outputValidateElement)
      _ <- checkIfStderrOutputExists(outputValidateElement)
    } yield outputValidateElement

    val graphSections: Vector[WorkflowGraphElement] = bodyElements.filterByType[WorkflowGraphElement]

    val metaSectionValidation: ErrorOr[Option[MetaSectionElement]] = validateSize(bodyElements.filterByType[MetaSectionElement], "meta", 1)
    val parameterMetaSectionValidation: ErrorOr[Option[ParameterMetaSectionElement]] = validateSize(bodyElements.filterByType[ParameterMetaSectionElement], "parameterMeta", 1)

    (inputsSectionValidation, outputsSectionValidation, metaSectionValidation, parameterMetaSectionValidation, intermediateValueDeclarationStdoutCheck, intermediateValueDeclarationStderrCheck) mapN {
      (validInputs, validOutputs, meta, parameterMeta, _, _) =>
      WorkflowDefinitionElement(name, validInputs, graphSections.toSet, validOutputs, meta, parameterMeta, sourceLocation)
    }
  }

  private def checkIfStdoutInputExist(inputSection: Option[InputsSectionElement]): ErrorOr[Option[String]] = {
    inputSection match {
      case Some(section) =>
        if (section.inputDeclarations.flatMap(_.expression).exists(_.isInstanceOf[StdoutElement.type])) {
          s"Workflow cannot have stdout expression in input section at workflow-level.".invalidNel
        } else None.validNel
      case None => None.validNel
    }
  }

  private def checkIfStdoutOutputExists(outputSection: Option[OutputsSectionElement]): ErrorOr[Option[String]] = {
    outputSection match {
      case Some(section) =>
        if (section.outputs.map(_.expression).exists(_.isInstanceOf[StdoutElement.type])) {
          s"Workflow cannot have stdout expression in output section at workflow-level.".invalidNel
        } else None.validNel
      case None => None.validNel
    }
  }

  private def checkStdoutIntermediates(intermediate: Vector[IntermediateValueDeclarationElement]): ErrorOr[Option[String]] = {
    if (intermediate.map(_.expression).exists(_.isInstanceOf[StdoutElement.type])) {
      s"Workflow cannot have stdout expression at intermediate declaration section at workflow-level.".invalidNel
    } else None.validNel
  }

  private def checkIfStderrInputExist(inputs: Option[InputsSectionElement]): ErrorOr[Option[String]] = {
    inputs match {
      case Some(section) =>
        if (section.inputDeclarations.flatMap(_.expression).exists(_.isInstanceOf[StderrElement.type])) {
          s"Workflow cannot have stderr expression in input section at workflow-level.".invalidNel
        } else None.validNel
      case None => None.validNel
    }
  }

  private def checkIfStderrOutputExists(outputs: Option[OutputsSectionElement]): ErrorOr[Option[String]] = {
    outputs match {
      case Some(section) =>
        if (section.outputs.map(_.expression).exists(_.isInstanceOf[StderrElement.type])) {
          s"Workflow cannot have stderr expression in output section at workflow-level.".invalidNel
        } else None.validNel
      case None => None.validNel
    }
  }

  private def checkStderrIntermediates(intermediates: Vector[IntermediateValueDeclarationElement]): ErrorOr[Option[String]] = {
    if (intermediates.map(_.expression).exists(_.isInstanceOf[StderrElement.type])) {
      s"Workflow cannot have stderr expression at intermediate declaration section at workflow-level.".invalidNel
    } else None.validNel
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
