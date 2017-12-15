package cwl

import shapeless._

object ArgumentToSortingKey extends Poly1 {
  implicit def script: Case.Aux[Expression, Option[Int]] = at[Expression] { _: Expression => Option.empty[Int] }

  implicit def clb: Case.Aux[CommandLineBinding, Option[Int]] = at[CommandLineBinding] { _.position }

  implicit def string: Case.Aux[String, Option[Int]] = at[String] { _ => Option.empty[Int] }
}
