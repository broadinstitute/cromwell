package cwl

import cats.syntax.either._
import cats.syntax.validated._
import cats.instances.list._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cwl.CwlVersion._
import cwl.ExpressionTool.{ExpressionToolInputParameter, ExpressionToolOutputParameter}
import shapeless.{Coproduct, Witness}
import wom.RuntimeAttributes
import wom.callable.{Callable, CallableExpressionTaskDefinition}
import wom.expression.{IoFunctionSet, ValueAsAnExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.types.{WomMapType, WomStringType, WomType}
import wom.values.{WomMap, WomObject, WomString, WomValue}

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
                           cwlVersion: Option[CwlVersion] = None,
                           `$namespaces`: Option[Map[String, String]] = None,
                           `$schemas`: Option[Array[String]] = None
                         ) extends Tool {

  def asCwl: Cwl = Coproduct[Cwl](this)

  def castWomValueAsMap(evaluatedExpression: WomValue): Checked[Map[String, WomValue]] = {
    evaluatedExpression match {
      case WomMap(WomMapType(WomStringType, _), map) =>  map.map{
        case (WomString(key), value) => key -> value
        case (key, _) => throw new RuntimeException(s"saw a non-string value $key in wom map $evaluatedExpression typed with string keys ")
      }.asRight
      case obj: WomObject => obj.values.asRight
      case _ => s"Could not cast value $evaluatedExpression to a Map[String, WomValue]".invalidNelCheck
    }
  }

  def buildTaskDefinition(taskName: String,
                          inputDefinitions: List[_ <: Callable.InputDefinition],
                          outputDefinitions: List[Callable.OutputDefinition],
                          runtimeAttributes: RuntimeAttributes,
                          requirementsAndHints: List[cwl.Requirement],
                          expressionLib: ExpressionLib): ErrorOr[CallableExpressionTaskDefinition] = {

    def evaluate(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, outputPorts: List[OutputPort]): Checked[Map[OutputPort, WomValue]] = {
      val womExpression = expression match {
        case StringOrExpression.String(str) => ValueAsAnExpression(WomString(str))
        case StringOrExpression.Expression(expr) => ECMAScriptWomExpression(expr, inputNames, expressionLib)
      }

      // If we expect a certain type for 
      def coerce(womValue: WomValue, womType: WomType): Checked[WomValue] = womType.coerceRawValue(womValue).toChecked

      /*
        * Look in the object (the result of the expression evaluation) for values to match output each port.
        * Ideally we'd want to find a value for each output port declared by the tool.
        * However since the spec is not clear on this and cwltool seems to ignore entirely the outputs declaration,
        * we'll do the same for now, meaning we'll assign a value to the output port iff the result of the expression has one for them,
        * otherwise we'll remove the port for the final map.
        * If there are fields in the expression's result that do not have a matching declared output port, they won't be added either,
        * because we'd have to create an output port for them now which would be useless anyway since nothing could point to it.
       */
      def mapPortsToValues(values: Map[String, WomValue]): Checked[Map[OutputPort, WomValue]] = {
        import cats.syntax.traverse._

        val coercedValues: List[ErrorOr[(OutputPort, WomValue)]] = for {
          outputPort <- outputPorts
          value <- values.get(outputPort.internalName)
          coerced = coerce(value, outputPort.womType).toValidated
        } yield coerced.map(outputPort -> _)

        coercedValues.sequence[ErrorOr, (OutputPort, WomValue)].map(_.toMap).toEither
      }

      for {
        // Evaluate the expression
        evaluatedExpression <- womExpression.evaluateValue(inputs, ioFunctionSet).toEither
        /*
          * We expect the result to be an object of the form:
          * {
          *   "outputName": value,
          *   ...
          * }
         */
        womMap <- castWomValueAsMap(evaluatedExpression)

        mappedValues <- mapPortsToValues(womMap)
      } yield mappedValues
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

object ExpressionTool {

  case class ExpressionToolInputParameter(id: String,
                                          label: Option[String] = None,
                                          secondaryFiles: Option[SecondaryFiles] = None,
                                          format: Option[InputParameterFormat] = None,
                                          streamable: Option[Boolean] = None,
                                          doc: Option[Doc] = None,
                                          inputBinding: Option[InputCommandLineBinding] = None,
                                          default: Option[CwlAny] = None,
                                          `type`: Option[MyriadInputType] = None) extends InputParameter

  case class ExpressionToolOutputParameter(id: String,
                                           label: Option[String] = None,
                                           secondaryFiles: Option[SecondaryFiles] = None,
                                           format: Option[OutputParameterFormat] = None,
                                           streamable: Option[Boolean] = None,
                                           doc: Option[Doc] = None,
                                           outputBinding: Option[CommandOutputBinding] = None,
                                           `type`: Option[MyriadOutputType] = None) extends OutputParameter
}
