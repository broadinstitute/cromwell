package wdl4s.cwl

import shapeless._
import wdl4s.wdl.command.StringCommandPart
import wdl4s.wom.CommandPart

object StringOrExpressionToCommandPart extends Poly1 {
  val EcmaScriptRegex = """\$\(([^)]*)\)""".r
  implicit def script = at[ECMAScript] { ex => CwlExpressionCommandPart(ex.value): CommandPart }

  implicit def fun = at[ECMAFunction] { f => CwlExpressionCommandPart(f.value): CommandPart }

  implicit def string = at[String] {
    case EcmaScriptRegex(expr) => CwlExpressionCommandPart(expr): CommandPart
    case part => StringCommandPart(part): CommandPart
  }
}
