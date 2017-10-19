package wom.callable

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wom.graph.{CallNode, GraphInputNode, LocalName, PortBasedGraphOutputNode}
import wom.callable.TaskDefinitionSpec._
import wom.types.{WdlIntegerType, WdlStringType}

class TaskDefinitionSpec extends FlatSpec with Matchers {

  // Checks that the graph generated from a task definition adds sufficient inputs and outputs, and is correctly linked.
  behavior of "TaskDefinition.graph"

  it should "represent an empty task as a one-node TaskCall graph" in {
    noInputsOrOutputsTask.graph match {
      case Valid(graph) =>
        graph.nodes.size should be(1)
      case Invalid(l) => fail(s"Failed to construct a one-node TaskCall graph: ${l.toList.mkString(", ")}")
    }
  }

  it should "create a graph input node for a one-input task definition" in {
    oneInputTask.graph match {
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
    oneOutputTask.graph match {
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
  
  it should "fail to build a graph with duplicates fqns" in {
    duplicateFqns.graph match {
      case Valid(_) => fail("The graph should be invalid")
      case Invalid(_) =>
    }
  }
}

object TaskDefinitionSpec {

  val noInputsOrOutputsTask = TaskDefinition(
    name = "foo",
    commandTemplate = Seq.empty,
    runtimeAttributes = null,
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List.empty,
    inputs = List.empty)

  val oneInputTask = TaskDefinition(
    name = "foo",
    commandTemplate = Seq.empty,
    runtimeAttributes = null,
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List.empty,
    inputs = List(Callable.RequiredInputDefinition(LocalName("bar"), WdlIntegerType)))

  val oneOutputTask = TaskDefinition(
    name = "foo",
    commandTemplate = Seq.empty,
    runtimeAttributes = null,
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List(Callable.OutputDefinition(LocalName("bar"), WdlStringType, null)),
    inputs = List.empty)
  
  val duplicateFqns = TaskDefinition(
    name = "foo",
    commandTemplate = Seq.empty,
    runtimeAttributes = null,
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List(Callable.OutputDefinition(LocalName("bar"), WdlStringType, null)),
    inputs = List(Callable.RequiredInputDefinition(LocalName("bar"), WdlStringType))
  )
}
