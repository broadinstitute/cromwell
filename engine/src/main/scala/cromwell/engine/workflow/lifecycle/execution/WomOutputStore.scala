package cromwell.engine.workflow.lifecycle.execution

import cromwell.core.CromwellGraphNode._
import cromwell.core.ExecutionIndex._
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.CollectorKey
import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.graph.GraphNode
import wdl4s.wom.graph.GraphNodePort.OutputPort

import scala.util.{Failure, Try}

object WomOutputStore {
  case class OutputEntry(name: String, wdlType: WdlType, wdlValue: Option[WdlValue])
  case class OutputCallKey(call: GraphNode, index: ExecutionIndex)
  def empty = OutputStore(Map.empty)
}

case class WomOutputStore(store: Map[OutputPort, WdlValue]) {

  override def toString = store.map { case (k, l) => s"$k -> $l" } mkString System.lineSeparator

  def add(values: Map[OutputPort, WdlValue]) = this.copy(store = store ++ values)

  def get(outputPort: OutputPort): Option[WdlValue] = store.get(outputPort)

//  def collectCall(call: TaskCallNode, scatter: ScatterNode, sortedShards: Seq[JobKey]) = Try {
//    val shardsOutputs = sortedShards map { e =>
//      fetchNodeOutputEntries(call, e.index) map {
//        case callOutputs: WdlCallOutputsObject => callOutputs.outputs
//        case _ => throw new RuntimeException("Call outputs should be a WdlCallOutputsObject")
//      } getOrElse(throw new RuntimeException(s"Could not retrieve output for shard ${e.scope} #${e.index}"))
//    }
//
//    call.callable.outputs map { taskOutput =>
//      val wdlValues = shardsOutputs.map(
//        _.getOrElse(taskOutput.unqualifiedName, throw new RuntimeException(s"Could not retrieve output ${taskOutput.unqualifiedName}")))
//      val arrayType = taskOutput.relativeWdlType(scatter).asInstanceOf[WdlArrayType]
//      val arrayOfValues = WdlArray(arrayType, wdlValues)
//      taskOutput.unqualifiedName -> JobOutput(arrayOfValues)
//    } toMap
//  }
//
//  def collectDeclaration(declaration: ExpressionNode, scatter: ScatterNode, sortedShards: Seq[JobKey]) = Try {
//    val shardsOutputs = sortedShards map { e =>
//      fetchNodeOutputEntries(declaration, e.index) getOrElse {
//        throw new RuntimeException(s"Could not retrieve output for shard ${e.scope} #${e.index}")
//      }
//    }
//    val arrayType = declaration.relativeWdlType(scatter).asInstanceOf[WdlArrayType]
//    Map(declaration.unqualifiedName -> JobOutput(WdlArray(arrayType, shardsOutputs)))
//  }

  /**
    * Try to generate output for a collector call, by collecting outputs for all of its shards.
    * It's fail-fast on shard output retrieval
    */
  def generateCollectorOutput(collector: CollectorKey,
                              shards: Iterable[JobKey]): Try[CallOutputs] = {
//    lazy val sortedShards = shards.toSeq sortBy { _.index.fromIndex }

    collector.scope match {
//      case call: TaskCallNode => collectCall(call, collector.scatter, sortedShards)
//      case declaration: ExpressionNode => collectDeclaration(declaration, collector.scatter, sortedShards)
      case other => Failure(new RuntimeException(s"Cannot retrieve outputs for ${other.fullyQualifiedName}"))
    }
  }
}
