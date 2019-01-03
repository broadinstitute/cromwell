package cwl.command

object ParentName {
  def empty: ParentName = ParentName(None)
  def apply(id: String): ParentName = ParentName(id.split("#").tail.headOption)
}

case class ParentName(value: Option[String]) {
  def stripParent(in: String) = value.map(v => in.stripPrefix(s"$v/")).getOrElse(in)
}
