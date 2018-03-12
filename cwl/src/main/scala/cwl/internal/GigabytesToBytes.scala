package cwl.internal

import shapeless.{Coproduct, Poly1}

object GigabytesToBytes extends Poly1 {
  val toMebibytesMultiplier = Math.pow(2, 20).toLong

  implicit def long = at[Long] {
    l =>
      val value = l * toMebibytesMultiplier
      Coproduct[cwl.ResourceRequirementType](value)
  }

  implicit def string = at[String] {
    s =>
      //TODO: Scale this by multiplier https://github.com/broadinstitute/cromwell/issues/3382
      Coproduct[cwl.ResourceRequirementType](s)
  }

  implicit def expression = at[cwl.Expression] {
    e =>
      //TODO: Scale this by multiplier https://github.com/broadinstitute/cromwell/issues/3382
      Coproduct[cwl.ResourceRequirementType](e)
  }
}

