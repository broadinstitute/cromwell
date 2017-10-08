
import wdl.values.WdlValue
import wom.callable.Callable.InputDefinition

package object wom {
  type WomEvaluatedCallInputs = Map[InputDefinition, WdlValue]
}
