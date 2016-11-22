package cromwell.engine.workflow.lifecycle.execution

import cromwell.core.ExecutionIndex._
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.OutputStore.{OutputCallKey, OutputEntry}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.CollectorKey
import wdl4s.types.{WdlArrayType, WdlType}
import wdl4s.util.TryUtil
import wdl4s.values.{WdlArray, WdlCallOutputsObject, WdlValue}
import wdl4s.{Call, Scope}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object OutputStore {
  case class OutputEntry(name: String, wdlType: WdlType, wdlValue: Option[WdlValue])
  case class OutputCallKey(call: Scope, index: ExecutionIndex)
  def empty = OutputStore(Map.empty)
}

case class OutputStore(store: Map[OutputCallKey, Traversable[OutputEntry]]) {
  def add(values: Map[OutputCallKey, Traversable[OutputEntry]]) = this.copy(store = store ++ values)

  def fetchCallOutputEntries(call: Call, index: ExecutionIndex): Try[WdlCallOutputsObject] = {
    def outputEntriesToMap(outputs: Traversable[OutputEntry]): Map[String, Try[WdlValue]] = {
      outputs map { output =>
        output.wdlValue match {
          case Some(wdlValue) => output.name -> Success(wdlValue)
          case None => output.name -> Failure(new RuntimeException(s"Could not retrieve output ${output.name} value"))
        }
      } toMap
    }

    store.get(OutputCallKey(call, index)) match {
      case Some(outputs) =>
        TryUtil.sequenceMap(outputEntriesToMap(outputs), s"Output fetching for call ${call.unqualifiedName}") map { outputsMap =>
          WdlCallOutputsObject(call, outputsMap)
        }
      case None => Failure(new RuntimeException(s"Could not find call ${call.unqualifiedName}"))
    }
  }

  /**
    * Try to generate output for a collector call, by collecting outputs for all of its shards.
    * It's fail-fast on shard output retrieval
    */
  def generateCollectorOutput(collector: CollectorKey,
                              shards: Iterable[CallKey]): Try[CallOutputs] = Try {
    val shardsOutputs = shards.toSeq sortBy { _.index.fromIndex } map { e =>
      fetchCallOutputEntries(e.scope, e.index) map {
        _.outputs
      } getOrElse(throw new RuntimeException(s"Could not retrieve output for shard ${e.scope} #${e.index}"))
    }
    collector.scope.outputs map { taskOutput =>
      val wdlValues = shardsOutputs.map(
        _.getOrElse(taskOutput.unqualifiedName, throw new RuntimeException(s"Could not retrieve output ${taskOutput.unqualifiedName}")))
      val arrayOfValues = new WdlArray(WdlArrayType(taskOutput.wdlType), wdlValues)
      taskOutput.unqualifiedName -> JobOutput(arrayOfValues)
    } toMap
  }
  
}
