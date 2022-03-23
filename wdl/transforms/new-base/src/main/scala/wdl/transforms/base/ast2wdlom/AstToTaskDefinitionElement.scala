package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
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
    val inputsSectionElement: ErrorOr[Option[InputsSectionElement]] = validateOneMax(bodyElements.collect { case e: InputsSectionElement => e }, "inputs")
    val declarations: Vector[IntermediateValueDeclarationElement] = bodyElements.collect { case e: IntermediateValueDeclarationElement => e }
    val outputsSectionElement: ErrorOr[Option[OutputsSectionElement]] = validateOneMax(bodyElements.collect { case e: OutputsSectionElement => e }, "outputs")
    val commandSectionElement: ErrorOr[CommandSectionElement] = validateExists(bodyElements.collect { case e: CommandSectionElement => e }, "command")
    val runtimeSectionElement: ErrorOr[Option[RuntimeAttributesSectionElement]] = validateOneMax(bodyElements.collect { case e: RuntimeAttributesSectionElement => e }, "runtime")

    val metaSectionElement: ErrorOr[Option[MetaSectionElement]] = validateOneMax(bodyElements.collect { case e: MetaSectionElement => e }, "meta")
    val parameterMetaSectionElement: ErrorOr[Option[ParameterMetaSectionElement]] = validateOneMax(bodyElements.collect { case e: ParameterMetaSectionElement => e }, "parameterMeta")

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
