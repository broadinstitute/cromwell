package wdl.draft3.transforms.wdlom2wom

import simulacrum.typeclass
import scala.language.implicitConversions

@typeclass
trait ValueGenerator[A] {
  def generatedValues(a: A): Set[String]
}

@typeclass
trait ValueConsumer[A] {
  def consumedValues(a: A): Set[ConsumedValue]
}

/**
  * Until we do the linking, we can't tell whether a consumed 'x.y' is a call output or a member access for 'y' on
  * a variable called 'x'.
  */
sealed trait ConsumedValue
final case class ConsumedSingleValue(name: String) extends ConsumedValue
final case class ConsumedLookupValue(name: String, firstLookup: String) extends ConsumedValue
