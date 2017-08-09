package wdltool

import wdl4s.wom.graph.GraphNodePort.{InputPort, OutputPort}
import wdl4s.wom.graph._

package object graph {

  private[graph] def dotSafe(s: String) = s""""${s.replaceAllLiterally("\"", "\\\"")}""""

  private[graph] implicit class GraphNodeGraphics(val graphNode: GraphNode) extends AnyVal {
    def graphFillColor = graphNode match {
      case _: GraphInputNode => "lightskyblue1"
      case _: GraphOutputNode => "palegreen"
      case _ => "white"
    }

    def graphName: String = dotSafe(graphNode match {
      case c: CallNode =>
        s"call ${c.name}"
      case s: ScatterNode =>
        s"scatter ${s.scatterVariableMapping.graphInputNode.name}"
      case gin: GraphInputNode =>
        s"${gin.womType.toWdlString} ${gin.name}"
      case gon: GraphOutputNode =>
        s"${gon.womType.toWdlString} ${gon.name}"
      case expr: ExpressionNode =>
        val inputNames = expr.instantiatedExpression.expression.inputs.mkString(", ")
        s"${expr.womType.toWdlString} ${expr.name} = f($inputNames)"
      case other =>
        throw new Exception(s"womgraph can't find a graphName for GraphNodes of type: ${other.getClass.getSimpleName}")
    })

    def graphId: String = dotSafe("NODE" + graphNode.hashCode().toString)
  }

  private[graph] implicit class GraphNodePortGraphics(val graphNodePort: GraphNodePort) extends AnyVal {
    def graphShape = graphNodePort match {
      case _: InputPort => "oval"
      case _: OutputPort => "hexagon"
    }

    def graphName: String = dotSafe(graphNodePort.name)
    def graphId: String = dotSafe("PORT" + graphNodePort.hashCode().toString)
  }
}
