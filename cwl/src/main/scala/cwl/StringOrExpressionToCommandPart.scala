package cwl

import shapeless._
import wdl.command.StringCommandPart

object StringOrExpressionToCommandPart extends Poly1 {
  implicit def script = at[Expression] { CwlExpressionCommandPart.apply _ }

  implicit def string = at[String] {
     StringCommandPart.apply _
  }
}
