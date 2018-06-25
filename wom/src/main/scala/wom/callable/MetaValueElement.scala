package wom.callable

sealed trait MetaValueElement

object MetaValueElement {
  // primitives: int, string, float, etc.
  final object MetaValueElementNull extends MetaValueElement
  final case class MetaValueElementBoolean(value: Boolean) extends MetaValueElement
  final case class MetaValueElementFloat(value: Double) extends MetaValueElement
  final case class MetaValueElementInteger(value: Int) extends MetaValueElement
  final case class MetaValueElementString(value: String) extends MetaValueElement

  // compounds: maps and arrays
  final case class MetaValueElementObject(value: Map[String, MetaValueElement]) extends MetaValueElement
  final case class MetaValueElementArray(value: Vector[MetaValueElement]) extends MetaValueElement
}

// A key-value pair, used in the meta sections
final case class MetaKvPair(key: String, value: MetaValueElement)
