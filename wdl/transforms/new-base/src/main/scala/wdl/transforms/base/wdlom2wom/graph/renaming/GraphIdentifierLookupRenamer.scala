package wdl.transforms.base.wdlom2wom.graph.renaming

import simulacrum.typeclass
import wdl.transforms.base.wdlom2wom.expression.renaming.IdentifierLookupRenamer
import wdl.model.draft3.elements.{ExpressionElement, WorkflowGraphElement}

@typeclass
trait GraphIdentifierLookupRenamer[A] {
  def renameIdentifiers(a: A, renamingMap: Map[String, String])
                       (implicit expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
                        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]): A
}
