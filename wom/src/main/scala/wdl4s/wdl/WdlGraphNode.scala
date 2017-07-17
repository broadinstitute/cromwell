package wdl4s.wdl

import wdl4s.wdl.AstTools.{EnhancedAstNode, VariableReference}

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

    (referencedNodes ++ closestScopedAncestor.toSeq).toSet
  }

  final lazy val downstream: Set[WdlGraphNode] = {
    for {
      node <- namespace.descendants.collect({
        case n: WdlGraphNode if n.fullyQualifiedName != fullyQualifiedName => n
      })
      if node.upstream.contains(this)
    } yield node
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
    variable <- expr.variableReferences
    scope <- parent
    node <- scope.resolveVariable(variable.terminal.sourceString)
  } yield node
}
