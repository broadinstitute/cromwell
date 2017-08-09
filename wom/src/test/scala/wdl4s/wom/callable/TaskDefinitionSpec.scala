package wdl4s.wom.callable

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.types.{WdlIntegerType, WdlStringType}
import wdl4s.wom.graph.{CallNode, GraphInputNode, PortBasedGraphOutputNode}

class TaskDefinitionSpec extends FlatSpec with Matchers {

  // Checks that the graph generated from a task definition adds sufficient inputs and outputs, and is correctly linked.
  behavior of "TaskDefinition.graph"

  it should "represent an empty task as a one-node TaskCall graph" in {
    val task = TaskDefinition(
      name = "foo",
      commandTemplate = Seq.empty,
      runtimeAttributes = null,
      meta = Map.empty,
      parameterMeta = Map.empty,
      outputs = Set.empty,
      inputs = Set.empty,
      declarations = List.empty)

    val graphValidation = task.graph

    graphValidation match {
      case Valid(graph) =>
        graph.nodes.size should be(1)
      case Invalid(l) => fail(s"Failed to construct a one-node TaskCall graph: ${l.toList.mkString(", ")}")
    }
  }

  it should "create a graph input node for a one-input task definition" in {
    val task = TaskDefinition(
      name = "foo",
      commandTemplate = Seq.empty,
      runtimeAttributes = null,
      meta = Map.empty,
      parameterMeta = Map.empty,
      outputs = Set.empty,
      inputs = Set(Callable.RequiredInputDefinition("bar", WdlIntegerType)),
      declarations = List.empty)

    val graphValidation = task.graph

    graphValidation match {
      case Valid(graph) =>
        graph.nodes.size should be(2)
        (graph.nodes.toList.find(_.isInstanceOf[GraphInputNode]), graph.nodes.toList.find(_.isInstanceOf[CallNode])) match {
          case (Some(inputNode), Some(callNode)) =>
            callNode.inputPorts.size should be(1)
            callNode.inputPorts.head.upstream.graphNode should be(inputNode)
          case other => fail(s"Oops: $other")
        }
      case Invalid(l) => fail(s"Failed to construct a one-input TaskCall graph: ${l.toList.mkString(", ")}")
    }
  }

  it should "create a graph output node for a one-output task definitions" in {
    val task = TaskDefinition(
      name = "foo",
      commandTemplate = Seq.empty,
      runtimeAttributes = null,
      meta = Map.empty,
      parameterMeta = Map.empty,
      outputs = Set(Callable.OutputDefinition("bar", WdlStringType, null)),
      inputs = Set(),
      declarations = List.empty)

    val graphValidation = task.graph

    graphValidation match {
      case Valid(graph) =>
        graph.nodes.size should be(2)
        (graph.nodes.toList.find(_.isInstanceOf[PortBasedGraphOutputNode]), graph.nodes.toList.find(_.isInstanceOf[CallNode])) match {
          case (Some(outputNode), Some(callNode)) =>
            callNode.outputPorts.size should be(1)
            outputNode.inputPorts.size should be(1)
            outputNode.inputPorts.head.upstream.graphNode should be(callNode)
          case other => fail(s"Oops: $other")
        }
      case Invalid(l) => fail(s"Failed to construct a one-input TaskCall graph: ${l.toList.mkString(", ")}")
    }
  }
}
