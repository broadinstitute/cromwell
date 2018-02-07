package wdl

import wdl.draft3.exception.OutputVariableLookupException
import wom.core.WorkflowSource
import wom.values.WomValue

import scala.util.{Failure, Try}

package object draft3 {
  type WorkflowJson = String
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WomValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type EvaluatedTaskInputs = Map[Declaration, WomValue]
  type ImportResolver = String => WorkflowSource
  type OutputResolver = (WdlGraphNode, Option[Int]) => Try[WomValue]

  val NoOutputResolver: OutputResolver = (node: WdlGraphNode, i: Option[Int]) => Failure(OutputVariableLookupException(node, i))
}
