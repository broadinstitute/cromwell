package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.either._
import common.collections.EnhancedCollections._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import wdl.model.draft3.elements._
import wom.SourceFileLocation

object AstToTaskDefinitionElement {

  def astToTaskDefinitionElement(implicit astNodeToTaskSectionElement: CheckedAtoB[GenericAstNode, TaskSectionElement]
                                ): CheckedAtoB[GenericAst, TaskDefinitionElement] = CheckedAtoB.fromErrorOr
  { a: GenericAst => s"read task definition${a.lineAndColumnString}" }
  { a =>

    val nameElementValidation: ErrorOr[String] = astNodeToString(a.getAttribute("name")).toValidated
    val sectionsValidation: ErrorOr[Vector[TaskSectionElement]] = a.getAttributeAsVector[TaskSectionElement]("sections").toValidated
    val sourceLocation : Option[SourceFileLocation] = a.getSourceLine.map(SourceFileLocation(_))

    (nameElementValidation, sectionsValidation) flatMapN { (nameElement, sections) =>
      combineElements(nameElement, sections, sourceLocation)
    }
  }

  def combineElements(nameElement: String,
                      bodyElements: Vector[TaskSectionElement],
                      sourceLocation: Option[SourceFileLocation]) = {
    val inputsSectionElement: ErrorOr[Option[InputsSectionElement]] = validateOneMax(bodyElements.filterByType[InputsSectionElement], "inputs")
    val declarations: Vector[IntermediateValueDeclarationElement] = bodyElements.filterByType[IntermediateValueDeclarationElement]
    val outputsSectionElement: ErrorOr[Option[OutputsSectionElement]] = validateOneMax(bodyElements.filterByType[OutputsSectionElement], "outputs")
    val commandSectionElement: ErrorOr[CommandSectionElement] = validateExists(bodyElements.filterByType[CommandSectionElement], "command")
    val runtimeSectionElement: ErrorOr[Option[RuntimeAttributesSectionElement]] = validateOneMax(bodyElements.filterByType[RuntimeAttributesSectionElement], "runtime")

    val metaSectionElement: ErrorOr[Option[MetaSectionElement]] = validateOneMax(bodyElements.filterByType[MetaSectionElement], "meta")
    val parameterMetaSectionElement: ErrorOr[Option[ParameterMetaSectionElement]] = validateOneMax(bodyElements.filterByType[ParameterMetaSectionElement], "parameterMeta")

    (inputsSectionElement, outputsSectionElement, commandSectionElement, runtimeSectionElement, metaSectionElement, parameterMetaSectionElement) mapN {
      (inputs, outputs, command, runtime, meta, parameterMeta) =>
        TaskDefinitionElement(nameElement, inputs, declarations, outputs, command, runtime, meta, parameterMeta, sourceLocation)
    }
  }

  private def validateExists[A](elements: Vector[A], sectionName: String): ErrorOr[A] = {
    val sectionValidation: ErrorOr[A] = if (elements.size != 1) {
      s"Task must have a $sectionName section.".invalidNel
    } else {
      elements.head.validNel
    }
    sectionValidation
  }

  private def validateOneMax[A](elements: Vector[A], sectionName: String): ErrorOr[Option[A]] = {
    val sectionValidation: ErrorOr[Option[A]] = if (elements.size > 1) {
      s"Task must have no more than one $sectionName sections, found ${elements.size}.".invalidNel
    } else {
      elements.headOption.validNel
    }
    sectionValidation
  }
}
