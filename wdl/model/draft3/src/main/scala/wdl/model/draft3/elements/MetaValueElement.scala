package wdl.model.draft3.elements

sealed trait MetaValueElement

object MetaValueElement {
  // primitives: int, string, float, etc.
  final object MNull extends MetaValueElement
  final case class MBoolean(value: Boolean) extends MetaValueElement
  final case class MFloat(value: Double) extends MetaValueElement
  final case class MInteger(value: Int) extends MetaValueElement
  final case class MString(value: String) extends MetaValueElement

  // compounds: maps and arrays
  final case class MObject(value: Map[String, MetaValueElement]) extends MetaValueElement
  final case class MArray(value: Vector[MetaValueElement]) extends MetaValueElement
}

// A key-value pair, used in the meta sections
final case class MetaKvPair(key: String, value: MetaValueElement)
