package wdl.draft2

import wdl.draft2.model.exception.OutputVariableLookupException
import wom.core.WorkflowSource
import wom.values.WomValue

import scala.util.{Failure, Try}

package object model {
  type WorkflowJson = String
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WomValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type EvaluatedTaskInputs = Map[Declaration, WomValue]
  type Draft2ImportResolver = String => WorkflowSource
  type OutputResolver = (WdlGraphNode, Option[Int]) => Try[WomValue]

  val NoOutputResolver: OutputResolver = (node: WdlGraphNode, i: Option[Int]) => Failure(OutputVariableLookupException(node, i))
}
