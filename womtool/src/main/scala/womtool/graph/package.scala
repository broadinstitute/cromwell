package womtool

import wom.graph.GraphNodePort.{InputPort, OutputPort}
import wom.graph._
import wom.graph.expression.{AnonymousExpressionNode, ExpressionNode}

package object graph {

  private[graph] def dotSafe(s: String) = s""""${s.replaceAllLiterally("\"", "\\\"")}""""

  private[graph] implicit class GraphNodeGraphics(val graphNode: GraphNode) extends AnyVal {
    def graphFillColor = graphNode match {
      case _: ConditionalNode | _: ScatterNode | _: WorkflowCallNode => "lightgray"
      case _: ExternalGraphInputNode => "lightskyblue1"
      case _: OuterGraphInputNode => "blueviolet"
      case _: PortBasedGraphOutputNode => "yellowgreen"
      case _: ExpressionBasedGraphOutputNode => "palegreen"
      case _ => "white"
    }

    def graphStyle = graphNode match {
      case _: AnonymousExpressionNode => "\"filled,dashed\""
      case o: OuterGraphInputNode if o.preserveScatterIndex => "\"dashed\""
      case _ => "\"filled,solid\""
    }

    def graphName: String = dotSafe(graphNode match {
      case c: CallNode => s"call ${c.fullyQualifiedName} (${c.localName})"
      case s: ScatterNode => s"scatter ${s.scatterCollectionExpressionNodes.head.identifier.localName.value} in"
      case _: ConditionalNode => "conditional"
      case gin: OptionalGraphInputNodeWithDefault => s"${gin.womType.stableName} ${gin.localName} = ..."
      case gin: GraphInputNode => s"${gin.womType.stableName} ${gin.localName}"
      case gon: GraphOutputNode => s"${gon.womType.stableName} ${gon.localName}"
      case expr: ExpressionNode =>
        val inputNames = expr.womExpression.inputs.mkString(", ")
        s"${expr.womType.stableName} ${expr.localName} = f($inputNames)"
      case other =>
        throw new Exception(s"womgraph can't find a graphName for GraphNodes of type: ${other.getClass.getSimpleName}")
    })

    def graphId: String = dotSafe("NODE" + graphObjectUniqueId(graphNode))
  }

  private[graph] implicit class GraphNodePortGraphics(val graphNodePort: GraphNodePort) extends AnyVal {
    def graphShape = graphNodePort match {
      case _: InputPort => "oval"
      case _: OutputPort => "hexagon"
    }

    def graphName: String = dotSafe(graphNodePort.womType.stableName + " " + graphNodePort.name)
    def graphId: String = dotSafe("PORT" + graphObjectUniqueId(graphNodePort))
  }

  /**
    * Should be good enough to provide a unique ID based on object reference.
    *
    * In some cases this breaks down (cf. https://stackoverflow.com/questions/10645494/can-i-assume-two-objects-with-the-same-system-identityhashcode-are-the-same#10645567),
    * but I think that should be rare enough to ignore (since womgraph is mostly a "help the developers" kind of a feature!)
    */
  private def graphObjectUniqueId(a: Any): Int = System.identityHashCode(a)
}
