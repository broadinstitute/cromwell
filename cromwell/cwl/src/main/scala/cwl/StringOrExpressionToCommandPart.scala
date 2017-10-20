package cwl

import cwl.command.StringCommandPart
import shapeless._

object StringOrExpressionToCommandPart extends Poly1 {
  implicit def script = at[Expression] { CwlExpressionCommandPart.apply }

  implicit def string = at[String] {
     StringCommandPart.apply
  }
}
