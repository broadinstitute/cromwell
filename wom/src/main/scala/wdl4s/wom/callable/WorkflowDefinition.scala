package wdl4s.wom.callable

import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wom.graph.Graph
import cats.syntax.validated._
import wdl4s.wom.expression.Expression
import wdl4s.wom.graph.GraphNode._

final case class WorkflowDefinition(name: String,
                                    innerGraph: Graph,
                                    meta: Map[String, String],
                                    parameterMeta: Map[String, String],
                                    declarations: List[(String, Expression)]) extends Callable {

  override lazy val toString = s"[Workflow $name]"
  override val graph: ErrorOr[Graph] = innerGraph.validNel

  override lazy val inputs: Set[_ <: Callable.InputDefinition] = innerGraph.nodes.inputDefinitions
}
