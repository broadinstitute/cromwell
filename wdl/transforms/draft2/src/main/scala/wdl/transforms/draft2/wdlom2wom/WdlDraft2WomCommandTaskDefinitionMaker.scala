package wdl.transforms.draft2.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft2.model.WdlTask
import wdl.draft2.model.{WdlTask, WdlWomExpression}
import wom.callable.Callable.{InputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.callable.{Callable, CallableTaskDefinition, CommandTaskDefinition}
import wom.graph.LocalName
import wom.transforms.WomCommandTaskDefinitionMaker
import wom.types.WomOptionalType

object WdlDraft2WomCommandTaskDefinitionMaker extends WomCommandTaskDefinitionMaker[WdlTask] {

  override def toWomTaskDefinition(wdlTask: WdlTask): ErrorOr[CommandTaskDefinition] = {

    val womInputs: List[Callable.InputDefinition] = (wdlTask.declarations collect {
      case d if d.expression.isEmpty && !d.womType.isInstanceOf[WomOptionalType] =>
        RequiredInputDefinition(LocalName(d.unqualifiedName), d.womType)
      case d if d.expression.isEmpty && d.womType.isInstanceOf[WomOptionalType] =>
        OptionalInputDefinition(LocalName(d.unqualifiedName), d.womType.asInstanceOf[WomOptionalType])
      case d if d.expression.nonEmpty =>
        InputDefinitionWithDefault(LocalName(d.unqualifiedName), d.womType, WdlWomExpression(d.expression.get, wdlTask))
    }).toList

    CallableTaskDefinition(
      name = wdlTask.fullyQualifiedName,
      commandTemplateBuilder = Function.const(wdlTask.commandTemplate.validNel),
      runtimeAttributes = wdlTask.runtimeAttributes.toWomRuntimeAttributes(wdlTask),
      meta = wdlTask.meta,
      parameterMeta = wdlTask.parameterMeta,
      outputs = wdlTask.outputs.map(_.womOutputDefinition).toList,
      inputs = womInputs,
      adHocFileCreation = Set.empty,
      environmentExpressions = Map.empty
    ).valid
  }
}
