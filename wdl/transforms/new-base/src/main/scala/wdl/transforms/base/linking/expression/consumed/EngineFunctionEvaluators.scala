package wdl.transforms.base.linking.expression.consumed

import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.{ExpressionValueConsumer, UnlinkedConsumedValueHook}
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._

object EngineFunctionEvaluators {

  implicit val stdoutElementValueConsumer: ExpressionValueConsumer[StdoutElement.type] = new ExpressionValueConsumer[ExpressionElement.StdoutElement.type] {
    override def expressionConsumedValueHooks(a: ExpressionElement.StdoutElement.type)
                                             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = Set.empty
  }

  implicit val stderrElementValueConsumer: ExpressionValueConsumer[StderrElement.type] = new ExpressionValueConsumer[ExpressionElement.StderrElement.type] {
    override def expressionConsumedValueHooks(a: ExpressionElement.StderrElement.type)
                                             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = Set.empty
  }

  implicit val readLinesValueConsumer: ExpressionValueConsumer[ReadLines] = forOneParamFunction
  implicit val readTsvValueConsumer: ExpressionValueConsumer[ReadTsv] = forOneParamFunction
  implicit val readMapValueConsumer: ExpressionValueConsumer[ReadMap] = forOneParamFunction
  implicit val readObjectValueConsumer: ExpressionValueConsumer[ReadObject] = forOneParamFunction
  implicit val readObjectsValueConsumer: ExpressionValueConsumer[ReadObjects] = forOneParamFunction
  implicit val readJsonValueConsumer: ExpressionValueConsumer[ReadJson] = forOneParamFunction
  implicit val readIntValueConsumer: ExpressionValueConsumer[ReadInt] = forOneParamFunction
  implicit val readStringValueConsumer: ExpressionValueConsumer[ReadString] = forOneParamFunction
  implicit val readFloatValueConsumer: ExpressionValueConsumer[ReadFloat] = forOneParamFunction
  implicit val readBooleanValueConsumer: ExpressionValueConsumer[ReadBoolean] = forOneParamFunction
  implicit val writeLinesValueConsumer: ExpressionValueConsumer[WriteLines] = forOneParamFunction
  implicit val writeTsvValueConsumer: ExpressionValueConsumer[WriteTsv] = forOneParamFunction
  implicit val writeMapValueConsumer: ExpressionValueConsumer[WriteMap] = forOneParamFunction
  implicit val writeObjectValueConsumer: ExpressionValueConsumer[WriteObject] = forOneParamFunction
  implicit val writeObjectsValueConsumer: ExpressionValueConsumer[WriteObjects] = forOneParamFunction
  implicit val writeJsonValueConsumer: ExpressionValueConsumer[WriteJson] = forOneParamFunction
  implicit val rangeValueConsumer: ExpressionValueConsumer[Range] = forOneParamFunction
  implicit val transposeValueConsumer: ExpressionValueConsumer[Transpose] = forOneParamFunction
  implicit val lengthValueConsumer: ExpressionValueConsumer[Length] = forOneParamFunction
  implicit val flattenValueConsumer: ExpressionValueConsumer[Flatten] = forOneParamFunction
  implicit val selectFirstValueConsumer: ExpressionValueConsumer[SelectFirst] = forOneParamFunction
  implicit val selectAllValueConsumer: ExpressionValueConsumer[SelectAll] = forOneParamFunction
  implicit val definedValueConsumer: ExpressionValueConsumer[Defined] = forOneParamFunction
  implicit val floorValueConsumer: ExpressionValueConsumer[Floor] = forOneParamFunction
  implicit val ceilValueConsumer: ExpressionValueConsumer[Ceil] = forOneParamFunction
  implicit val roundValueConsumer: ExpressionValueConsumer[Round] = forOneParamFunction
  implicit val globValueConsumer: ExpressionValueConsumer[Glob] = forOneParamFunction

  implicit val sizeValueConsumer: ExpressionValueConsumer[Size] = forOneOrTwoParamFunction
  implicit val basenameValueConsumer: ExpressionValueConsumer[Basename] = forOneOrTwoParamFunction

  implicit val zipValueConsumer: ExpressionValueConsumer[Zip] = forTwoParamFunction
  implicit val crossValueConsumer: ExpressionValueConsumer[Cross] = forTwoParamFunction
  implicit val prefixValueConsumer: ExpressionValueConsumer[Prefix] = forTwoParamFunction
  implicit val joinValueConsumer: ExpressionValueConsumer[Sep] = forTwoParamFunction

  implicit val subFunctionValueConsumer: ExpressionValueConsumer[Sub] = new ExpressionValueConsumer[Sub] {
    override def expressionConsumedValueHooks(a: Sub)
                                             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = {
      a.input.expressionConsumedValueHooks ++ a.pattern.expressionConsumedValueHooks ++ a.replace.expressionConsumedValueHooks
    }
  }

  private def forOneParamFunction[A <: OneParamFunctionCallElement]: ExpressionValueConsumer[A] = new ExpressionValueConsumer[A] {
    override def expressionConsumedValueHooks(a: A)
                                             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = a.param.expressionConsumedValueHooks
  }

  private def forOneOrTwoParamFunction[A <: OneOrTwoParamFunctionCallElement]: ExpressionValueConsumer[A] = new ExpressionValueConsumer[A] {
    override def expressionConsumedValueHooks(a: A)
                                             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = {
      a.firstParam.expressionConsumedValueHooks ++
        a.secondParam.toSet.flatMap { secondParam: ExpressionElement => secondParam.expressionConsumedValueHooks }
    }
  }

  private def forTwoParamFunction[A <: TwoParamFunctionCallElement]: ExpressionValueConsumer[A] = new ExpressionValueConsumer[A] {
    override def expressionConsumedValueHooks(a: A)
                                             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): Set[UnlinkedConsumedValueHook] = {
      a.arg1.expressionConsumedValueHooks ++ a.arg2.expressionConsumedValueHooks
    }
  }
}
