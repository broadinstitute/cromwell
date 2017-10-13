
import wdl.exception.OutputVariableLookupException
import wom.values.WdlValue
import wom.core._

import scala.util.{Failure, Try}

package object wdl {
  type WorkflowJson = String
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WdlValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type EvaluatedTaskInputs = Map[Declaration, WdlValue]
  type ImportResolver = String => WorkflowSource
  type OutputResolver = (WdlGraphNode, Option[Int]) => Try[WdlValue]

  val NoOutputResolver: OutputResolver = (node: WdlGraphNode, i: Option[Int]) => Failure(OutputVariableLookupException(node, i))
}
