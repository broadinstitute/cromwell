import wdl4s.values.WdlValue

import scala.util.{Failure, Try}

package object wdl4s {
  type WdlSource = String
  type WdlJson = String
  type WorkflowRawInputs = Map[FullyQualifiedName, Any]
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WdlValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type EvaluatedTaskInputs = Map[Declaration, WdlValue]
  type ImportResolver = String => WdlSource
  type OutputResolver = (GraphNode, Option[Int]) => Try[WdlValue]
  
  class OutputVariableLookupException(node: GraphNode, index: Option[Int]) extends VariableLookupException(s"Could not find outputs for call ${node.fullyQualifiedName} at index $index")
  val NoOutputResolver: OutputResolver = (node: GraphNode, i: Option[Int]) => Failure(new OutputVariableLookupException(node, i))

  trait TsvSerializable {
    def tsvSerialize: Try[String]
  }

}
