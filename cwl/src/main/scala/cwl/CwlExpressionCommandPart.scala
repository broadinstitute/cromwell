package cwl

import shapeless.{Inl, Inr}
import wdl.values.{WdlSingleFile, WdlString, WdlValue}
import wom.CommandPart
import wom.expression.IoFunctionSet
import wom.graph.LocalName


case class CwlExpressionCommandPart(expr: Expression) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WdlValue],
                           functions: IoFunctionSet,
                           valueMapper: (WdlValue) => WdlValue): String = {

    val pc = ParameterContext.Empty.withInputs(inputsMap.map({ case (LocalName(localName), value) => localName -> value }), functions)

    val wdlValue: WdlValue = expr.fold(EvaluateExpression).apply(pc)

    wdlValue.valueString
  }
}

case class CwlArgumentCommandPart(argument: CommandLineBinding) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WdlValue],
                           functions: IoFunctionSet,
                           valueMapper: (WdlValue) => WdlValue) = {

    val pc = ParameterContext.Empty.withInputs(inputsMap.map({
      case (LocalName(localName), WdlSingleFile(path)) => localName -> WdlString(path)
      case (LocalName(localName), value) => localName -> value
    }), functions)

    val wdlValue: WdlValue= argument match {
      case CommandLineBinding(_, _, _, _, _, Some(Inl(expression: Expression)), Some(false)) =>
        expression.fold(EvaluateExpression).apply(pc)
      case CommandLineBinding(_, _, _, _, _, Some(Inr(Inl(string: String))), Some(false)) =>
        WdlString(string)
      case _ => ???
    }

    wdlValue.valueString
  }
}

