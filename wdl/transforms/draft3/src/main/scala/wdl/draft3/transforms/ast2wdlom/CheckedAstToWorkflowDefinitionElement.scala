package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import common.validation.Checked._
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{InputsSectionElement, OutputsSectionElement, WorkflowDefinitionElement}

object CheckedAstToWorkflowDefinitionElement {

  def convert(a: Ast): ErrorOr[WorkflowDefinitionElement] = {

    val nameElementValidation: ErrorOr[String] = CheckedAtoB[AstNode, String].run(a.getAttribute("name")).toValidated

    val inputsAndOutputs: Checked[(Vector[InputsSectionElement], Vector[OutputsSectionElement])] =
      a.getAttributeAsVectors[InputsSectionElement, OutputsSectionElement]("body", "Inputs", "Outputs")

    val inputsValidation: ErrorOr[Option[InputsSectionElement]] = (inputsAndOutputs flatMap { case (i, _) =>
      if (i.size <= 1) { i.headOption.validNelCheck } else s"Workflow must have one or zero output sections, but got ${i.size} instead.".invalidNelCheck
    }).toValidated

    val outputsValidation: ErrorOr[Vector[OutputsSectionElement]] = (inputsAndOutputs map { _._2 }).toValidated


    (nameElementValidation, inputsValidation, outputsValidation) mapN { (name, inputs, outputs) => WorkflowDefinitionElement(name, inputs, outputs) }
  }
}
