package cwl

import cats.Applicative
import cats.data.Validated.{Invalid, Valid}
import cats.data.{Kleisli, NonEmptyList, OptionT}
import cats.data.Kleisli._
import cats.data.Validated._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.Checked
import common.validation.ErrorOr._
import cwl.CommandLineTool._
import cwl.command.ParentName
import shapeless.Coproduct
import wom.types.WomStringType
import wom.values.WomEvaluatedCallInputs

object CommandPartSortingAlgorithm {
  type HasShellCommandRequirement = Boolean

  type Inputs = (SchemaDefRequirement, ExpressionLib, WomEvaluatedCallInputs, HasShellCommandRequirement)

  type CommandPartFunc[A] = Kleisli[ErrorOr, (SchemaDefRequirement, ExpressionLib, WomEvaluatedCallInputs, HasShellCommandRequirement), A]

  def inputBindingsCommandParts(inputs: Array[CommandInputParameter]): CommandPartFunc[List[SortKeyAndCommandPart]] =
    inputs.toList.flatTraverse[CommandPartFunc, SortKeyAndCommandPart](inputBindingsCommandPart)

  def inputBindingsCommandPart(inputParameter: CommandInputParameter): CommandPartFunc[List[SortKeyAndCommandPart]] =
    Kleisli{ case ((schemaDefRequirement, expressionLib, inputValues, hasShellCommandRequirement)) =>
      val parsedName = FullyQualifiedName(inputParameter.id)(ParentName.empty).id

      val womType = inputParameter.`type`.map(_.fold(MyriadInputTypeToWomType).apply(schemaDefRequirement)).getOrElse(WomStringType)


      val defaultValue = inputParameter.default.map(_.fold(InputParameter.DefaultToWomValuePoly).apply(womType))

      inputValues
        .collectFirst({ case (inputDefinition, womValue) if inputDefinition.name == parsedName => womValue.validNel })
        .orElse(defaultValue) match {
        case Some(Valid(value)) =>
          // See http://www.commonwl.org/v1.0/CommandLineTool.html#Input_binding
          lazy val initialKey = CommandBindingSortingKey.empty
            .append(inputParameter.inputBinding, Coproduct[StringOrInt](parsedName))

          inputParameter.`type`.toList.
            flatMap(
              _.fold(MyriadInputTypeToSortedCommandParts).
                apply(inputParameter.inputBinding, value, initialKey.asNewKey, hasShellCommandRequirement, expressionLib, schemaDefRequirement)).validNel
        case Some(Invalid(errors)) => Invalid(errors)
        case None => s"Could not find an input value for input $parsedName in ${inputValues.prettyString}".invalidNel
      }
    }

}
