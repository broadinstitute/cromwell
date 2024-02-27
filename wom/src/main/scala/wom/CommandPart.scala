package wom

import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WomValue

trait CommandPart {
  def instantiate(inputsMap: Map[LocalName, WomValue],
                  functions: IoFunctionSet,
                  valueMapper: WomValue => WomValue
  ): ErrorOr[List[InstantiatedCommand]]
}
