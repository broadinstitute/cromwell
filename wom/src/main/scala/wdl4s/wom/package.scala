package wdl4s

import wdl4s.wdl.values.WdlValue
import wdl4s.wom.callable.Callable.InputDefinition

package object wom {
  type WomEvaluatedCallInputs = Map[InputDefinition, WdlValue]
}
