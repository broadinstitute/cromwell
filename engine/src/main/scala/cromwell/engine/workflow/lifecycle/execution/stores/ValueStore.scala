package cromwell.engine.workflow.lifecycle.execution.stores

import cats.syntax.validated._
import cats.syntax.option._
import cromwell.core.ExecutionIndex._
import cromwell.engine.workflow.lifecycle.execution.keys.{ConditionalCollectorKey, ScatterCollectorKey}
import common.collections.Table
import common.validation.ErrorOr.ErrorOr
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore.ValueKey
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.values.WomArray.WomArrayLike
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

  override def toString: String = {
    val values = store.valuesTriplet.map {
      case (node, index, value) => s"$node:${index.fromIndex} -> $value"
    }

    s"""
      |ValueStore {
      |  ${values.mkString("," + System.lineSeparator + "  ")}
      |}
    """.stripMargin
  }

  final def add(values: Map[ValueKey, WomValue]): ValueStore = {
    this.copy(store = store.addAll(values.map({ case (key, value) => (key.port, key.index, value) })))
  }

  final def get(outputKey: ValueKey): Option[WomValue] = store.getValue(outputKey.port, outputKey.index)

  final def get(outputPort: OutputPort, index: ExecutionIndex): Option[WomValue] = store.getValue(outputPort, index)

  /**
    * Collect all shards values for a node
    */
  final def collectShards(collector: ScatterCollectorKey): ErrorOr[Map[ValueKey, WomArray]] = {
    /*
     * collector.scatter.outputMapping is all the output ports for the scatter (ScatterGathererPort)
     * There is one for each PortBasedGraphOutputNode in the scatter inner graph
     */
    val scatterGatherPort = collector.scatterGatherPort

    // The output port which this scatterGatherPort is referring to (could be a call node output port or a declaration node output port for example)
    val sourcePort = scatterGatherPort.outputToGather.source

    def goodResult(values: Seq[WomValue]) = Map(ValueKey(scatterGatherPort, None) -> WomArray(scatterGatherPort.womType, values)).validNel

    if (collector.scatterWidth == 0) {
      goodResult(Seq.empty)
    } else {
      // Try to find the values for this port in the store. Remember, the value store is a Table[OutputPort, ExecutionIndex, WomValue],
      // so by "get"ting the output port we get all the shards for this output port, which is exactly what we want here (we're collecting the shards)
      store.rowOptional(sourcePort) match {
        case Some(shards) =>
          val collectedValue = shards.toList.sortBy(_._1).map(_._2)

          // Verify that we have all the expected shards.
          if (collectedValue.size == collector.scatterWidth) {
            // If the sizes match, create an Array from the values and assign it to the ScatterGatherPort
            goodResult(collectedValue)
          } else {
            // If we don't find enough shards, this collector was found "runnable" when it shouldn't have
            s"Some shards are missing from the value store for node ${collector.node.fullyQualifiedName}, expected ${collector.scatterWidth} shards but only got ${collectedValue.size}: ${collectedValue.mkString(", ")}".invalidNel
          }
        case None =>
          // If we don't find any shards, this collector was found "runnable" when it shouldn't have
          s"Scatter collector cannot find any values for ${sourcePort.identifier.fullyQualifiedName.value} in value store: $this".invalidNel
      }
    }
  }

  final def collectConditional(collector: ConditionalCollectorKey): ErrorOr[Map[ValueKey, WomOptionalValue]] = {
    val conditionalPort = collector.conditionalOutputPort
    val sourcePort = conditionalPort.outputToExpose.source
    store.getValue(sourcePort, collector.index) match {
      case Some(womValue) => Map(ValueKey(conditionalPort, collector.index) -> WomOptionalValue(womValue).flattenOptional).validNel
      case None => s"Conditional collector cannot find a value for output port ${sourcePort.identifier.fullyQualifiedName.value} in value store: $this".invalidNel
    }
  }

  def resolve(index: Option[Int])(outputPort: OutputPort): ErrorOr[WomValue] = {

    // Look up the value of the scatter variable's expression Node at the appropriate index
    def forScatterVariable(svn: ScatterVariableNode): ErrorOr[WomValue] = {
      val p: OutputPort = svn.linkToOuterGraph
      get(p, None) match {
        // Try to find the element at "jobIndex" in the array value stored for the outputPort, any other case is a failure
        case Some(womValue: WomArrayLike) =>
          index match {
            case Some(jobIndex) =>
              womValue.asArray.value.lift(jobIndex)
                .toValidNel(s"Shard index $jobIndex exceeds scatter array length: ${womValue.asArray.value.size}")
            case None => s"Unsharded execution key references a scatter variable: ${p.identifier.fullyQualifiedName}".invalidNel
          }
        case Some(other) => s"Value for scatter collection ${p.identifier.fullyQualifiedName} is not an array: ${other.womType.toDisplayString}".invalidNel
        case None => s"Can't find a value for scatter collection ${p.identifier.fullyQualifiedName} (looking for index $index)".invalidNel
      }
    }

    // Just look it up in the store
    def forGraphNodePort(p: OutputPort, index: ExecutionIndex) = get(p, index) match {
      case Some(value) => value.validNel
      case None =>
        s"Can't find a ValueStore value for $p at index $index in $this".invalidNel
    }

    def findValueStorePort(p: OutputPort, index: ExecutionIndex): ErrorOr[WomValue] = {
      p.graphNode match {
        case svn: ScatterVariableNode => forScatterVariable(svn)
        case ogin: OuterGraphInputNode if ogin.preserveScatterIndex => findValueStorePort(ogin.linkToOuterGraph, index)
        case ogin: OuterGraphInputNode => findValueStorePort(ogin.linkToOuterGraph, None)
        case _: GraphInputNode => forGraphNodePort(p, None) // Must be a workflow input, which never have indices
        case _ => forGraphNodePort(p, index)
      }
    }

    findValueStorePort(outputPort, index)
  }
}
