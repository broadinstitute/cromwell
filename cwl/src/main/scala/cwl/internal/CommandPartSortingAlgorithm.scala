package cwl.internal

import cats.data.Kleisli._
import cats.data.ReaderT
import cats.data.Validated._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import cwl.CommandLineTool._
import cwl.command.ParentName
import cwl.{ArgumentCommandLineBinding, ArgumentToCommandPart, CommandLineTool, CommandPartExpression, FullyQualifiedName, InputParameter, MyriadInputTypeToSortedCommandParts, MyriadInputTypeToWomType}
import shapeless.Coproduct
import wom.types.WomStringType

object CommandPartSortingAlgorithm {
  def argumentCommandParts(arguments: Option[Array[CommandLineTool.Argument]]): CommandPartExpression[List[SortKeyAndCommandPart]] =
      // arguments is an Option[Array[Argument]], the toList.flatten gives a List[Argument]
    arguments.toList.flatten
      // zip the index because we need it in the sorting key
      .zipWithIndex.flatTraverse(argumentToCommandPart.tupled)

  def argumentToCommandPart: (Argument, Int) => CommandPartExpression[List[SortKeyAndCommandPart]] = (argument, index) => ReaderT {
    case ((requirementsAndHints, expressionLib, _)) =>
      val part = argument.fold(ArgumentToCommandPart).apply(requirementsAndHints.hasShellCommandRequirement, expressionLib)
      // Get the position from the binding if there is one
      val position = argument.select[ArgumentCommandLineBinding].flatMap(_.position)
        .map(Coproduct[StringOrInt](_)).getOrElse(CommandLineTool.DefaultPosition)

      // The key consists of the position followed by the index
      val sortingKey = CommandBindingSortingKey(List(position, Coproduct[StringOrInt](index)))

      List(SortKeyAndCommandPart(sortingKey, part)).validNel

  }

  def inputBindingsCommandParts(inputs: Array[CommandInputParameter]): CommandPartExpression[List[SortKeyAndCommandPart]] =
    inputs.toList.flatTraverse(inputBindingsCommandPart)

  def inputBindingsCommandPart(inputParameter: CommandInputParameter): CommandPartExpression[List[SortKeyAndCommandPart]] =
    ReaderT{ case ((hintsAndRequirements, expressionLib, inputValues)) =>
      val parsedName = FullyQualifiedName(inputParameter.id)(ParentName.empty).id

      val womType = inputParameter.`type`.map(_.fold(MyriadInputTypeToWomType).apply(hintsAndRequirements.schemaDefRequirement)).getOrElse(WomStringType)

      val defaultValue = inputParameter.default.map(_.fold(InputParameter.DefaultToWomValuePoly).apply(womType))

      inputValues
        .collectFirst({ case (inputDefinition, womValue) if inputDefinition.name == parsedName => womValue.validNel })
        .orElse(defaultValue) match {
        case Some(Valid(value)) =>
          // See http://www.commonwl.org/v1.0/CommandLineTool.html#Input_binding
          lazy val initialKey = CommandBindingSortingKey.empty
            .append(inputParameter.inputBinding, Coproduct[StringOrInt](parsedName))

          inputParameter.`type`.toList.
            flatMap{
              _.fold(MyriadInputTypeToSortedCommandParts).
                apply(
                  inputParameter.inputBinding,
                  value,
                  initialKey.asNewKey,
                  hintsAndRequirements.hasShellCommandRequirement,
                  expressionLib,
                  hintsAndRequirements.schemaDefRequirement)
            }.validNel
        case Some(Invalid(errors)) => Invalid(errors)
        case None => s"Could not find an input value for input $parsedName in ${inputValues.prettyString}".invalidNel
      }
    }

}
