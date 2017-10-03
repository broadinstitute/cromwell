package cwl

import scala.Function._
import shapeless._
import wdl.command.StringCommandPart

object ArgumentToCommandPart extends Poly1 {
  implicit def script = at[Expression] { const(StringCommandPart(null)) }

  implicit def clb = at[CommandLineBinding] {
      //TODO: This option.get will not hold up under scrutiny
      _.valueFrom.map(_.fold(StringOrExpressionToCommandPart)).get

      //TODO: Shell Quote = false?
  }

  implicit def string = at[String] {
    StringCommandPart.apply
  }
}
