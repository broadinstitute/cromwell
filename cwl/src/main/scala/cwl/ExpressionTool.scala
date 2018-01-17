package cwl

import cats.syntax.validated._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import cwl.CwlVersion._
import shapeless.{:+:, CNil, Coproduct, Witness}
import wom.RuntimeAttributes
import wom.callable.{Callable, CallableExpressionTaskDefinition}
import wom.expression.{IoFunctionSet, ValueAsAnExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.values.{WomString, WomValue}

case class ExpressionTool(
                           inputs: Array[ExpressionToolInputParameter] = Array.empty,
                           outputs: Array[ExpressionToolOutputParameter] = Array.empty,
                           `class`: Witness.`"ExpressionTool"`.T,
                           expression: StringOrExpression,
                           id: String,
                           requirements: Option[Array[Requirement]] = None,
                           hints: Option[Array[Hint]] = None,
                           label: Option[String] = None,
                           doc: Option[String] = None,
                           cwlVersion: Option[CwlVersion] = None) extends Tool {

  def asCwl: Cwl = Coproduct[Cwl](this)

  def buildTaskDefinition(taskName: String,
                          inputDefinitions: List[_ <: Callable.InputDefinition],
                          outputDefinitions: List[Callable.OutputDefinition],
                          runtimeAttributes: RuntimeAttributes,
                          requirementsAndHints: List[cwl.Requirement],
                          expressionLib: ExpressionLib): ErrorOr[CallableExpressionTaskDefinition] = {

    def evaluate(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, outputPorts: List[OutputPort]): Checked[Map[OutputPort, WomValue]] = {
      val womExpression = expression match {
        case StringOrExpression.String(str) => ValueAsAnExpression(WomString(str))
        case StringOrExpression.Expression(expr) => JavascriptExpression(expr, inputNames, expressionLib)
      }
      
      ExpressionToolEvaluation.evaluate(inputs, ioFunctionSet, outputPorts, womExpression)
    }

    CallableExpressionTaskDefinition(
      taskName,
      evaluate,
      runtimeAttributes,
      Map.empty,
      Map.empty,
      outputDefinitions,
      inputDefinitions,
      // TODO: This doesn't work in all cases and it feels clunky anyway - find a way to sort that out
      prefixSeparator = "#"
    ).validNel
  }
}

case class ExpressionToolOutputParameter(id: String,
                                         label: Option[String] = None,
                                         secondaryFiles: Option[SecondaryFiles] = None,
                                         format: Option[StringOrExpression] = None,
                                         streamable: Option[Boolean] = None,
                                         doc: Option[String :+: Array[String] :+: CNil] = None,
                                         outputBinding: Option[CommandOutputBinding] = None,
                                         `type`: Option[MyriadOutputType]) extends OutputParameter

case class ExpressionToolInputParameter(id: String,
                                        label: Option[String] = None,
                                        secondaryFiles: Option[SecondaryFiles] = None,
                                        format: Option[Expression :+: String :+: Array[String] :+: CNil] = None,
                                        streamable: Option[Boolean] = None,
                                        doc: Option[String :+: Array[String] :+: CNil] = None,
                                        inputBinding: Option[InputCommandLineBinding] = None,
                                        default: Option[CwlAny] = None,
                                        `type`: Option[MyriadInputType] = None) extends InputParameter
