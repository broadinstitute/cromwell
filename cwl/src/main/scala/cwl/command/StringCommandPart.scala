package cwl.command

import wom.CommandPart
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WomValue

case class StringCommandPart(literal: String) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WomValue], functions: IoFunctionSet, valueMapper: (WomValue) => WomValue) = literal
}
