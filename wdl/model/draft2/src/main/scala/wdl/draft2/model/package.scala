package wdl.draft2

import wdl.draft2.model.exception.OutputVariableLookupException
import wom.ResolvedImportRecord
import wom.core.WorkflowSource
import wom.values.WomValue

import scala.util.{Failure, Try}

case class Draft2ResolvedImportBundle(source: WorkflowSource, resolvedImportRecord: ResolvedImportRecord)

package object model {
  type WorkflowJson = String
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WomValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type EvaluatedTaskInputs = Map[Declaration, WomValue]
  type Draft2ImportResolver = String => Draft2ResolvedImportBundle
  type OutputResolver = (WdlGraphNode, Option[Int]) => Try[WomValue]

  val NoOutputResolver: OutputResolver = (node: WdlGraphNode, i: Option[Int]) => Failure(OutputVariableLookupException(node, i))
}
