package wdl

import cats.data.Validated.Valid
import cats.syntax.validated._
import common.collections.EnhancedCollections._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl.AstTools.{EnhancedAstNode, VariableReference}
import wom.graph.CallNode.CallNodeAndNewNodes
import wom.graph.GraphNode.GeneratedNodeAndNewNodes
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExposedExpressionNode


sealed trait WdlGraphNode extends Scope {

  /**
    * The set of all graph nodes which are (transitively) upstream from this one.
    */
  final lazy val upstreamAncestry: Set[WdlGraphNode] = WdlGraphNode.calculateUpstreamAncestry(Set.empty, this)

  def referencedNodes: Iterable[WdlGraphNode]

  /**
    * The set of all graph nodes which are a single step upstream from this one.
    */
  final lazy val upstream: Set[WdlGraphNode] = {
    // If we are inside the scope of another graph node (i.e. via document element ancestry), then
    // that is also upstream of us.
    val closestScopedAncestor = ancestry.collectFirst({
      case ancestor: WdlGraphNode => ancestor
    })

    // We want:
    // - Nodes that this node references
    // - The immediate ancestor
    // - Anything that our children sre downstream from
    // But because our children's upstream might also include these (which we don't want), filter out:
    // - This
    // - Any other WdlGraphNode descendants of this
    (referencedNodes ++ closestScopedAncestor.toSeq ++ childGraphNodes.flatMap(_.upstream)).toSet - this -- descendants.filterByType[WdlGraphNode]
  }

  final lazy val downstream: Set[WdlGraphNode] = {
    for {
      node <- namespace.descendants.collect({
        case n: WdlGraphNode if n.fullyQualifiedName != fullyQualifiedName => n
      })
      if node.upstream.contains(this)
    } yield node
  }

  def isUpstreamFrom(other: WdlGraphNode): Boolean = {
    other.upstreamAncestry.contains(this) || (other.childGraphNodes exists isUpstreamFrom)
  }
}

object WdlGraphNode {
  private def calculateUpstreamAncestry(currentSet: Set[WdlGraphNode], graphNode: WdlGraphNode): Set[WdlGraphNode] = {
    val setWithUpstream = currentSet ++ graphNode.upstream
    val updatesNeeded = graphNode.upstream -- currentSet
    updatesNeeded.foldLeft(setWithUpstream)(calculateUpstreamAncestry)
  }
}

trait WdlGraphNodeWithUpstreamReferences extends WdlGraphNode {
  def upstreamReferences: Iterable[VariableReference]

  // If we have variable reference to other graph nodes, then they are upstream from us.
  override final def referencedNodes = for {
      variable <- upstreamReferences
      node <- resolveVariable(variable.terminal.sourceString)
      if node.fullyQualifiedNameWithIndexScopes != fullyQualifiedNameWithIndexScopes
    } yield node
}

trait WdlGraphNodeWithInputs extends WdlGraphNode {
  def inputMappings: Map[String, WdlExpression]

  override final def referencedNodes = for {
    expr <- inputMappings.values
    variable <- expr.variableReferences(this)
    scope <- parent
    node <- scope.resolveVariable(variable.terminal.sourceString)
  } yield node
}
