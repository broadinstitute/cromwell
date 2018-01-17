package cwl

import cats.syntax.either._
import common.Checked
import common.validation.Validation._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.types.{WomObjectType, WomType}
import wom.values.{WomObjectLike, WomValue}

object ExpressionToolEvaluation {
  def evaluate(inputs: Map[String, WomValue],
               ioFunctionSet: IoFunctionSet,
               outputPorts: List[OutputPort],
               womExpression: WomExpression): Checked[Map[OutputPort, WomValue]] = {

    // Since we don't validate anything, if something goes wrong just remove it from the list of outputs, hence the toOption
    def coerce(womValue: WomValue, womType: WomType): Option[WomValue] = womType.coerceRawValue(womValue).toOption

    def mapPortsToValues(values: Map[String, WomValue]): Map[OutputPort, WomValue] = {
      outputPorts.flatMap({ outputPort =>
        values.get(outputPort.name)
          .flatMap(coerce(_, outputPort.womType))
          .map(outputPort -> _)
      }).toMap
    }

    for {
      evaluatedExpression <- womExpression.evaluateValue(inputs, ioFunctionSet).toEither
      asObject <- WomObjectType.coerceRawValue(evaluatedExpression).toChecked
      values = asObject match {
        case womObject: WomObjectLike => womObject.values
        case _ => throw new Exception("This cannot happen otherwise the coercion above would have failed")
      }
    } yield mapPortsToValues(values)
  }
}
