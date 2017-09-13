package cromwell.engine.workflow.lifecycle.execution

import cromwell.core.CromwellGraphNode._
import cromwell.core.ExecutionIndex._
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.OutputStore.OutputKey
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.CollectorKey
import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.graph.GraphNodePort.OutputPort
import wdl4s.wom.graph._

import scala.util.{Failure, Try}

object OutputStore {
  case class OutputEntry(name: String, wdlType: WdlType, wdlValue: Option[WdlValue])
  case class OutputKey(port: OutputPort, index: ExecutionIndex)
  def empty = OutputStore(Map.empty)
}

case class OutputStore(store: Map[OutputKey, WdlValue]) {

  override def toString = store.map { case (k, l) => s"$k -> ${l.valueString}" } mkString System.lineSeparator

  def add(values: Map[OutputKey, WdlValue]) = this.copy(store = store ++ values)

  def get(outputKey: OutputKey): Option[WdlValue] = store.get(outputKey)
  
  def get(outputPort: OutputPort, index: ExecutionIndex): Option[WdlValue] = store.get(OutputKey(outputPort, index))

  // TODO WOM: re-implement the OutputStore for WOM
  def generateCollectorOutput(collector: CollectorKey,
                              shards: Iterable[JobKey]): Try[CallOutputs] = {
    //lazy val sortedShards = shards.toSeq sortBy { _.index.fromIndex }
    
    collector.scope match {
      case _: CallNode => Failure(new NotImplementedError("Scatters are not working in WOMland yet"))
      case _: ExpressionNode => Failure(new NotImplementedError("Scatters are not working in WOMland yet"))
      case other => Failure(new RuntimeException(s"Cannot retrieve outputs for ${other.fullyQualifiedName}")) 
    }
  }
}
