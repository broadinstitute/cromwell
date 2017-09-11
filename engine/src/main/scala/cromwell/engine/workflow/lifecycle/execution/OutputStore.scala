package cromwell.engine.workflow.lifecycle.execution

import cromwell.core.CromwellGraphNode._
import cromwell.core.ExecutionIndex._
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.OutputStore.{OutputCallKey, OutputEntry}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.CollectorKey
import lenthall.util.TryUtil
import wdl4s.wdl._
import wdl4s.wdl.types.{WdlArrayType, WdlType}
import wdl4s.wdl.values.{WdlArray, WdlCallOutputsObject, WdlValue}
import wdl4s.wom.graph.{ExpressionNode, GraphNode, ScatterNode, TaskCallNode}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object OutputStore {
  case class OutputEntry(name: String, wdlType: WdlType, wdlValue: Option[WdlValue])
  case class OutputCallKey(call: GraphNode, index: ExecutionIndex)
  def empty = OutputStore(Map.empty)
}

case class OutputStore(store: Map[OutputCallKey, List[OutputEntry]]) {

  // TODO WOM: temporary conversion until WomOutputStore becomes THE OutputStore
  def toWomOutputStore: WomOutputStore = {
    val nodeMap = store.flatMap {
      case (key, values) =>
        key.call.outputPorts map { outputPort => outputPort -> values.find(_.name == outputPort.name).get.wdlValue.get }
    }
    WomOutputStore(nodeMap)
  }
  
  override def toString = store.map { case (k, l) => s"$k -> ${l.mkString(" ")}" } mkString System.lineSeparator

  def add(values: Map[OutputCallKey, List[OutputEntry]]) = this.copy(store = store ++ values)

  def fetchNodeOutputEntries(node: GraphNode, index: ExecutionIndex): Try[WdlValue] = {
    def outputEntriesToMap(outputs: List[OutputEntry]): Map[String, Try[WdlValue]] = {
      outputs map { output =>
        output.wdlValue match {
          case Some(wdlValue) => output.name -> Success(wdlValue)
          case None => output.name -> Failure(new RuntimeException(s"Could not retrieve output ${output.name} value"))
        }
      } toMap
    }

    def callOutputs(call: WdlCall, outputs: List[OutputEntry]) = {
      TryUtil.sequenceMap(outputEntriesToMap(outputs), s"Output fetching for call ${node.unqualifiedName}") map { outputsMap =>
        WdlCallOutputsObject(call, outputsMap)
      }
    }
    
    def declarationOutputs(declaration: Declaration, outputs: List[OutputEntry]) = {
      outputs match {
        case OutputEntry(_, _, Some(value)) :: Nil => Success(value)
        case _ => Failure(new RuntimeException(s"Could not find value for declaration ${declaration.fullyQualifiedName}"))
      }
    }
    
    store.get(OutputCallKey(node, index)) match {
      case Some(outputs) =>
        node match {
          case call: WdlCall => callOutputs(call, outputs)
          case declaration: Declaration => declarationOutputs(declaration, outputs)
          case other =>  Failure(new RuntimeException(s"Only Calls and Declarations are allowed in the OutputStore, found ${other.getClass.getSimpleName}"))
        }
      case None => Failure(new RuntimeException(s"Could not find scope ${node.unqualifiedName}"))
    }
  }
  
  def collectCall(call: TaskCallNode, scatter: ScatterNode, sortedShards: Seq[JobKey]) = Try {
    val shardsOutputs = sortedShards map { e =>
      fetchNodeOutputEntries(call, e.index) map {
        case callOutputs: WdlCallOutputsObject => callOutputs.outputs
        case _ => throw new RuntimeException("Call outputs should be a WdlCallOutputsObject")
      } getOrElse(throw new RuntimeException(s"Could not retrieve output for shard ${e.scope} #${e.index}"))
    }
    
    call.callable.outputs map { taskOutput =>
      val wdlValues = shardsOutputs.map(
        _.getOrElse(taskOutput.unqualifiedName, throw new RuntimeException(s"Could not retrieve output ${taskOutput.unqualifiedName}")))
      val arrayType = taskOutput.relativeWdlType(scatter).asInstanceOf[WdlArrayType]
      val arrayOfValues = WdlArray(arrayType, wdlValues)
      taskOutput.unqualifiedName -> JobOutput(arrayOfValues)
    } toMap
  }
  
  def collectDeclaration(declaration: ExpressionNode, scatter: ScatterNode, sortedShards: Seq[JobKey]) = Try {
    val shardsOutputs = sortedShards map { e =>
      fetchNodeOutputEntries(declaration, e.index) getOrElse {
        throw new RuntimeException(s"Could not retrieve output for shard ${e.scope} #${e.index}")
      }
    }
    val arrayType = declaration.relativeWdlType(scatter).asInstanceOf[WdlArrayType]
    Map(declaration.unqualifiedName -> JobOutput(WdlArray(arrayType, shardsOutputs)))
  }

  /**
    * Try to generate output for a collector call, by collecting outputs for all of its shards.
    * It's fail-fast on shard output retrieval
    */
  def generateCollectorOutput(collector: CollectorKey,
                              shards: Iterable[JobKey]): Try[CallOutputs] = {
    lazy val sortedShards = shards.toSeq sortBy { _.index.fromIndex }
    
    collector.scope match {
      case call: TaskCallNode => collectCall(call, collector.scatter, sortedShards)
      case declaration: ExpressionNode => collectDeclaration(declaration, collector.scatter, sortedShards)
      case other => Failure(new RuntimeException(s"Cannot retrieve outputs for ${other.fullyQualifiedName}")) 
    }
  }
}
