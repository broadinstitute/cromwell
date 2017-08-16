package wdl4s.wdl.command
import wdl4s.wdl.{Declaration, EvaluatedTaskInputs}
import wdl4s.wdl.expression.WdlFunctions
import wdl4s.wdl.values.WdlValue

case class CwlExpressionCommandPart(expr: String) extends CommandPart {
  override def instantiate(
                            declarations: Seq[Declaration],
                            inputsMap: EvaluatedTaskInputs,
                            functions: WdlFunctions[WdlValue],
                            valueMapper: (WdlValue) => WdlValue): String = ???
}

