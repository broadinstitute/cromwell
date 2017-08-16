package wdl4s.cwl

import shapeless._
import wdl4s.wdl.command.{CommandPart, CwlExpressionCommandPart, StringCommandPart}

object StringOrExpressionToCommandPart extends Poly1 {

  implicit def expression = at[ECMAScriptExpression] { ex => CwlExpressionCommandPart(ex.value): CommandPart }

  implicit def string = at[String] {
    StringCommandPart(_): CommandPart
  }
}
