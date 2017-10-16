package cromwell.engine.workflow.lifecycle.execution

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import cromwell.core.ExecutionIndex._
import cromwell.engine.workflow.lifecycle.execution.ValueStore.ValueKey
import cromwell.engine.workflow.lifecycle.execution.keys.{ConditionalCollectorKey, ScatterCollectorKey}
import lenthall.collections.Table
import lenthall.validation.ErrorOr.ErrorOr
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.values.{WomArray, WomOptionalValue, WomValue}

object ValueStore {
  def initialize(knownValues: Map[OutputPort, WomValue]): ValueStore = {
    val wdlValues = knownValues map { case (port, value) => (port, None, value) }
    ValueStore(Table.fill(wdlValues))
  }

  case class ValueKey(port: OutputPort, index: ExecutionIndex)
  def empty = ValueStore(Table.empty)
}

case class ValueStore(store: Table[OutputPort, ExecutionIndex, WomValue]) {

  override def toString: String = store.valuesTriplet.map {
    case (node, index, value) => s"${node.identifier.fullyQualifiedName}:${index.fromIndex} -> ${value.valueString}"
  } mkString System.lineSeparator

  final def add(values: Map[ValueKey, WomValue]): ValueStore = {
    this.copy(store = store.addAll(values.map({ case (key, value) => (key.port, key.index, value) })))
  }

  final def get(outputKey: ValueKey): Option[WomValue] = store.getValue(outputKey.port, outputKey.index)

  final def get(outputPort: OutputPort, index: ExecutionIndex): Option[WomValue] = store.getValue(outputPort, index)

  /**
    * Collect all shards values for a node
    */
  final def collectShards(collector: ScatterCollectorKey): ErrorOr[Map[ValueKey, WomArray]] = {
    lazy val outputValues: ErrorOr[Map[ValueKey, WomArray]] =
    /*
     * collector.scatter.outputMapping is all the output ports for the scatter (ScatterGathererPort)
     * There is one for each PortBasedGraphOutputNode in the scatter inner graph
     * (FIXME: ExpressionBasedGraphOutputNodes should probably be in there too for sub workflows)
    */
      collector.scatterGatherPorts.toList
        .traverse[ErrorOr, (ValueKey, WomArray)]({ scatterGatherPort =>
        // The output port which this scatterGatherPort is referring to (could be a call node output port or a declaration node output port for example)
        val sourcePort = scatterGatherPort.outputToGather.source

        // Try to find the values for this port in the store. Remember, the value store is a Table[OutputPort, ExecutionIndex, WomValue],
        // so by "get"ting the output port we get all the shards for this output port, which is exactly what we want here (we're collecting the shards)
        store.rowOptional(sourcePort) match {
          case Some(shards) =>
            val collectedValue = shards.toList.sortBy(_._1).map(_._2)

            // Verify that we have all the expected shards.
            if (collectedValue.size == collector.scatterWidth) {
              // If the sizes match, create an Array from the values and assign it to the ScatterGatherPort
              (ValueKey(scatterGatherPort, None) -> WomArray(scatterGatherPort.womType, collectedValue)).validNel
            } else {
              //If not something went really wrong and this collector was found "runnable" when it shouldn't have
              s"Some shards are missing from the value store, expected ${collector.scatterWidth} shards but only got ${collectedValue.size}: ${collectedValue.mkString(", ")}".invalidNel
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

  final def collectConditional(collector: ConditionalCollectorKey, isInBypassed: Boolean): ErrorOr[Map[ValueKey, WomOptionalValue]] = {
    val conditionalPort = collector.conditionalOutputPort
    val sourcePort = conditionalPort.outputToExpose.source
    store.getValue(sourcePort, collector.index) match {
      case Some(womValue: WomOptionalValue) if isInBypassed && womValue.value.isEmpty => Map(ValueKey(conditionalPort, collector.index) -> womValue).validNel
      case Some(womValue) if !isInBypassed  => Map(ValueKey(conditionalPort, collector.index) -> WomOptionalValue(womValue)).validNel
      case Some(_)  => s"Found a non empty value in a bypassed conditional node ${sourcePort.identifier.fullyQualifiedName.value}".invalidNel
      case None => s"Cannot find a value for output port ${sourcePort.identifier.fullyQualifiedName.value}".invalidNel
    }
  }
}
