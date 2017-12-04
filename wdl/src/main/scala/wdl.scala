
import wdl.exception.OutputVariableLookupException
import wom.values.WomValue
import wom.core._

import scala.util.{Failure, Try}

package object wdl {
  type WorkflowJson = String
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WomValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type EvaluatedTaskInputs = Map[Declaration, WomValue]
  type ImportResolver = String => WorkflowSource
  type OutputResolver = (WdlGraphNode, Option[Int]) => Try[WomValue]

  val NoOutputResolver: OutputResolver = (node: WdlGraphNode, i: Option[Int]) => Failure(OutputVariableLookupException(node, i))
}
