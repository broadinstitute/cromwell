package wom.callable

import cats.syntax.validated._
import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wom.graph.{CallNode, GraphInputNode, LocalName, PortBasedGraphOutputNode}
import wom.callable.CommandTaskDefinitionSpec._
import wom.types.{WomIntegerType, WomStringType}

class CommandTaskDefinitionSpec extends FlatSpec with Matchers {

  // Checks that the graph generated from a task definition adds sufficient inputs and outputs, and is correctly linked.
  behavior of "TaskDefinition.graph"

  it should "represent an empty task as a one-node TaskCall graph" in {
    executableNoInputsOrOutputsTask match {
      case Valid(task) =>
        task.graph.nodes.size should be(1)
      case Invalid(l) => fail(s"Failed to construct a one-node TaskCall graph: ${l.toList.mkString(", ")}")
    }
  }

  it should "create a graph input node for a one-input task definition" in {
    executableOneInputTask match {
      case Valid(task) =>
        task.graph.nodes.size should be(2)
        (task.graph.nodes.toList.find(_.isInstanceOf[GraphInputNode]), task.graph.nodes.toList.find(_.isInstanceOf[CallNode])) match {
          case (Some(inputNode), Some(callNode)) =>
            callNode.inputPorts.size should be(1)
            callNode.inputPorts.head.upstream.graphNode should be(inputNode)
          case other => fail(s"Oops: $other")
        }
      case Invalid(l) => fail(s"Failed to construct a one-input TaskCall graph: ${l.toList.mkString(", ")}")
    }
  }

  it should "create a graph output node for a one-output task definitions" in {
    executableOneOutputTask match {
      case Valid(task) =>
        val graph = task.graph
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
    executableDuplicateFqns match {
      case Valid(_) => fail("The graph should be invalid")
      case Invalid(e) =>
        e.size should be(1)
        e.head should be("Two or more nodes have the same FullyQualifiedName: foo.bar")
    }
  }
}

object CommandTaskDefinitionSpec {

  val noInputsOrOutputsTask = CallableTaskDefinition(
    name = "foo",
    commandTemplateBuilder = Function.const(Seq.empty.validNel),
    runtimeAttributes = null,
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List.empty,
    inputs = List.empty,
    adHocFileCreation = Set.empty,
    environmentExpressions = Map.empty,
    sourceLocation = None)
  val executableNoInputsOrOutputsTask = noInputsOrOutputsTask.toExecutable

  val oneInputTask = CallableTaskDefinition(
    name = "foo",
    commandTemplateBuilder = Function.const(Seq.empty.validNel),
    runtimeAttributes = null,
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List.empty,
    inputs = List(Callable.RequiredInputDefinition(LocalName("bar"), WomIntegerType)),
    adHocFileCreation = Set.empty,
    environmentExpressions = Map.empty,
    sourceLocation = None)
  val executableOneInputTask = oneInputTask.toExecutable

  val oneOutputTask = CallableTaskDefinition(
    name = "foo",
    commandTemplateBuilder = Function.const(Seq.empty.validNel),
    runtimeAttributes = null,
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List(Callable.OutputDefinition(LocalName("bar"), WomStringType, null)),
    inputs = List.empty,
    adHocFileCreation = Set.empty,
    environmentExpressions = Map.empty,
    sourceLocation = None)
  val executableOneOutputTask = oneOutputTask.toExecutable

  val duplicateFqns = CallableTaskDefinition(
    name = "foo",
    commandTemplateBuilder = Function.const(Seq.empty.validNel),
    runtimeAttributes = null,
    meta = Map.empty,
    parameterMeta = Map.empty,
    outputs = List(Callable.OutputDefinition(LocalName("bar"), WomStringType, null)),
    inputs = List(Callable.RequiredInputDefinition(LocalName("bar"), WomStringType), Callable.RequiredInputDefinition(LocalName("bar"), WomStringType)),
    adHocFileCreation = Set.empty,
    environmentExpressions = Map.empty,
    sourceLocation = None
  )
  val executableDuplicateFqns = duplicateFqns.toExecutable
}
