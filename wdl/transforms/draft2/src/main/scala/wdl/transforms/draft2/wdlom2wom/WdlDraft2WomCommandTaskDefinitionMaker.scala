package wdl.transforms.draft2.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft2.model.{AstTools, WdlTask, WdlWomExpression}
import wdl.draft2.parser.WdlParser.Terminal
import wom.SourceFileLocation
import wom.callable.Callable.{OverridableInputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.callable.{Callable, CallableTaskDefinition, CommandTaskDefinition, MetaValueElement}
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
        OverridableInputDefinitionWithDefault(LocalName(d.unqualifiedName), d.womType, WdlWomExpression(d.expression.get, wdlTask))
    }).toList

    // Figure out the start line of the workflow in the source file
    val t: Option[Terminal] = AstTools.findTerminals(wdlTask.ast).headOption

    // Draft-2 only support string values. It does not support composite values, or
    // anything
    def stringifyMetaValues(meta: Map[String, String]): Map[String, MetaValueElement] = {
      meta map {
        case (key, value) =>
          key -> MetaValueElement.MetaValueElementString(value)
      }
    }

    CallableTaskDefinition(
      name = wdlTask.fullyQualifiedName,
      commandTemplateBuilder = Function.const(wdlTask.commandTemplate.validNel),
      runtimeAttributes = wdlTask.runtimeAttributes.toWomRuntimeAttributes(wdlTask),
      meta = stringifyMetaValues(wdlTask.meta),
      parameterMeta = stringifyMetaValues(wdlTask.parameterMeta),
      outputs = wdlTask.outputs.map(_.womOutputDefinition).toList,
      inputs = womInputs,
      adHocFileCreation = Set.empty,
      environmentExpressions = Map.empty,
      sourceLocation = t.map(x => SourceFileLocation(x.getLine))
    ).valid
  }
}
