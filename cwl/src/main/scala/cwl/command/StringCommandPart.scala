package cwl.command

import wom.CommandPart
import wom.expression.IoFunctionSet
import wom.graph.LocalName
import wom.values.WdlValue

case class StringCommandPart(literal: String) extends CommandPart {
  override def instantiate(inputsMap: Map[LocalName, WdlValue], functions: IoFunctionSet, valueMapper: (WdlValue) => WdlValue) = literal
}
