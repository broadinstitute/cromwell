package cromwell.engine.workflow.lifecycle.execution

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import cromwell.core.ExecutionIndex._
import cromwell.engine.workflow.lifecycle.execution.ValueStore.OutputKey
import cromwell.engine.workflow.lifecycle.execution.keys.{ConditionalCollectorKey, ScatterCollectorKey}
import lenthall.collections.EnhancedCollections._
import lenthall.collections.Table
import lenthall.collections.Table.Table
import lenthall.validation.ErrorOr.ErrorOr
import wdl.types.WdlType
import wdl.values.{WdlArray, WdlOptionalValue, WdlValue}
import wom.graph.GraphNodePort.OutputPort
import wom.graph._

object ValueStore {
  def initialize(knownValues: Map[OutputPort, WdlValue]): ValueStore = {
    val wdlValues = knownValues map { case (port, value) => (port, None, value) }
    ValueStore(Table.fill(wdlValues))
  }

  case class OutputEntry(name: String, wdlType: WdlType, wdlValue: Option[WdlValue])
  case class OutputKey(port: OutputPort, index: ExecutionIndex)
  def empty = ValueStore(Map.empty)
}

case class ValueStore(store: Table[OutputPort, ExecutionIndex, WdlValue]) {

  override def toString: String = store.valuesTriplet.map {
    case (node, index, value) => s"${node.identifier.fullyQualifiedName}:${index.fromIndex} -> ${value.valueString}"
  } mkString System.lineSeparator

  final def add(values: Map[OutputKey, WdlValue]): ValueStore = {
    this.copy(store = store.addAll(values.map({ case (key, value) => (key.port, key.index, value) })))
  }

  final def get(outputKey: OutputKey): Option[WdlValue] = store.get(outputKey.port, outputKey.index)

  final def get(outputPort: OutputPort, index: ExecutionIndex): Option[WdlValue] = store.get(outputPort, index)

  /**
    * Collect all shards values for a node
    */
  final def collectShards(collector: ScatterCollectorKey): ErrorOr[Map[OutputKey, WdlArray]] = {
    lazy val outputValues: ErrorOr[Map[OutputKey, WdlArray]] =
    /*
     * collector.scatter.outputMapping is all the output ports for the scatter (ScatterGathererPort)
     * There is one for each PortBasedGraphOutputNode in the scatter inner graph 
     * (FIXME: ExpressionBasedGraphOutputNodes should probably be in there too for sub workflows)
    */
      collector.scatterGatherPorts.toList
        .traverse[ErrorOr, (OutputKey, WdlArray)]({ scatterGatherPort =>
        // The output port which this scatterGatherPort is referring to (could be a call node output port or a declaration node output port for example)
        val sourcePort = scatterGatherPort.outputToGather.source

        // Try to find the values for this port in the store. Remember, the value store is a Table[OutputPort, ExecutionIndex, WdlValue],
        // so by "get"ting the output port we get all the shards for this output port, which is exactly what we want here (we're collecting the shards)
        store.get(sourcePort) match {
          case Some(shards) =>
            val collectedValue = shards.toList.sortBy(_._1).map(_._2)

            // Verify that we have all the expected shards.
            if (collectedValue.size == collector.scatterWidth) {
              // If the sizes match, create an Array from the values and assign it to the ScatterGatherPort
              (OutputKey(scatterGatherPort, None) -> WdlArray(scatterGatherPort.womType, collectedValue)).validNel
            } else {
              //If not something went really wrong and this collector was found "runnable" when it shouldn't have
              "Some shards are missing from the value store".invalidNel
            }
          case None =>
            //If we don't something went really wrong and this collector was found "runnable" when it shouldn't have
            s"Cannot find a value for ${sourcePort.identifier.fullyQualifiedName.value} in the value store".invalidNel
        }
      })
        // If the traversing is successful convert the list of pairs to a Map
        .map(_.toMap)

    collector.node match {
      case _: CallNode | _: ExpressionNode | _: ConditionalNode => outputValues
      case other => s"Cannot retrieve outputs for ${other.fullyQualifiedName}".invalidNel
    }
  }

  final def collectConditional(collector: ConditionalCollectorKey, isInBypassed: Boolean): ErrorOr[Map[OutputKey, WdlOptionalValue]] = {
    val conditionalPort = collector.conditionalOutputPort
    val sourcePort = conditionalPort.outputToExpose.source
    store.get(sourcePort, collector.index) match {
      case Some(wdlValue: WdlOptionalValue) if isInBypassed && wdlValue.value.isEmpty => Map(OutputKey(conditionalPort, collector.index) -> wdlValue).validNel
      case Some(wdlValue) if !isInBypassed  => Map(OutputKey(conditionalPort, collector.index) -> WdlOptionalValue(wdlValue)).validNel
      case Some(_)  => s"Found a non empty value in a bypassed conditional node ${sourcePort.identifier.fullyQualifiedName.value}".invalidNel
      case None => s"Cannot find a value for output port ${sourcePort.identifier.fullyQualifiedName.value}".invalidNel
    }
  }
}
