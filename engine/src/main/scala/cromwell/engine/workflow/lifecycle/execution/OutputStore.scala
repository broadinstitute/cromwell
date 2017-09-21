package cromwell.engine.workflow.lifecycle.execution

import cromwell.core.CromwellGraphNode._
import cromwell.core.ExecutionIndex._
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.OutputStore.OutputKey
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.CollectorKey
import wdl.types.WdlType
import wdl.values.WdlValue
import wom.executable.Executable.ResolvedExecutableInputs
import wom.graph.GraphNodePort.OutputPort
import wom.graph._

import scala.util.{Failure, Try}

object OutputStore {
  def initialize(knownValues: ResolvedExecutableInputs): OutputStore = {
    val wdlValues = knownValues flatMap {
      // Known wdl values can be added to the output store
      case (port, resolvedValue) => resolvedValue.select[WdlValue] map { OutputKey(port, None) -> _ }
    }
    OutputStore(wdlValues)
  }

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
    
    collector.node match {
      case _: CallNode => Failure(new NotImplementedError("Scatters are not working in WOMland yet"))
      case _: ExpressionNode => Failure(new NotImplementedError("Scatters are not working in WOMland yet"))
      case other => Failure(new RuntimeException(s"Cannot retrieve outputs for ${other.fullyQualifiedName}")) 
    }
  }
}
