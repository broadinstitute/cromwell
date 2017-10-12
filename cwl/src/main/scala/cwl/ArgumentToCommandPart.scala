package cwl

import shapeless._
import wdl.command.StringCommandPart

object ArgumentToCommandPart extends Poly1 {
  implicit def script = at[Expression] { CwlExpressionCommandPart.apply }

  implicit def clb = at[CommandLineBinding] { CwlArgumentCommandPart.apply }

  implicit def string = at[String] { StringCommandPart.apply }
}
