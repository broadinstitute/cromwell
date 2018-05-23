package wdl.draft3.transforms.wdlom2wom.graph.renaming

import simulacrum.typeclass
import wdl.draft3.transforms.wdlom2wom.expression.renaming.IdentifierLookupRenamer
import wdl.model.draft3.elements.{ExpressionElement, WorkflowGraphElement}

import scala.language.implicitConversions

@typeclass
trait GraphIdentifierLookupRenamer[A] {
  def renameIdentifiers(a: A, renamingMap: Map[String, String])
                       (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
                        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]): A
}
