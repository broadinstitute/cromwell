package wdl.draft3.transforms.wdlom2wom

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.either._
import cats.syntax.traverse._
import cats.instances.list._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{CommandPartElement, FileBodyElement, TaskDefinitionElement}
import wom.{CommandPart, RuntimeAttributes}
import wom.callable.{Callable, CallableTaskDefinition, TaskDefinition}


object TaskDefinitionElementToWomTaskDefinition {
  def convert(a: FileBodyElement): ErrorOr[TaskDefinition] = a match {
    case a: TaskDefinitionElement =>
      // TODO inputs, outputs, runtime attributes
      val validInputs: ErrorOr[List[Callable.InputDefinition]] = List.empty.valid
      val validOutputs: ErrorOr[List[Callable.OutputDefinition]] = List.empty.valid
      val validRuntimeAttributes: ErrorOr[RuntimeAttributes] = RuntimeAttributes(Map.empty).valid

      implicit val commandPartConverter: CheckedAtoB[CommandPartElement, CommandPart] = commandPartElementToWomCommandPart

      val validCommand: ErrorOr[Seq[CommandPart]] = {
        a.commandSection.parts.toList.traverse[ErrorOr, CommandPart] { parts =>
          commandPartConverter.run(parts).toValidated
        }.map(_.toSeq)
      }

      (validInputs, validOutputs, validRuntimeAttributes, validCommand) mapN { (inputs, outputs, runtime, command) =>
        CallableTaskDefinition(a.name, Function.const(command.validNel), runtime, Map.empty, Map.empty, outputs, inputs, Set.empty, Map.empty)
      }

  }
}
