package wom.graph

/**
  * Identifies a node in its local context.
  * Used in during graph construction to link nodes together and validate the graph.
  */
case class LocalName(value: String) {
  def combineToLocalName(other: String) = LocalName(s"$value.$other")
  def combineToFullyQualifiedName(other: String): FullyQualifiedName = FullyQualifiedName(s"$value.$other")
  def combineToFullyQualifiedName(other: LocalName): FullyQualifiedName = combineToFullyQualifiedName(other.value)
}

/**
  * Identifies a node in its graph.
  * It must be unique for all CallNodes, GraphInputNodes and GraphOutputNodes within a graph.
  * Otherwise the Graph won't validate.
  * This *can* be the same as local name.
  * It is not required by WOM strictly speaking but rather useful to implementations
  * for serializing or reporting.
  */
case class FullyQualifiedName(value: String) {
  def combine(other: String) = FullyQualifiedName(s"$value.$other")
  def prefixWith(prefix: String) = FullyQualifiedName(s"$prefix.$value")
}

object WomIdentifier {
  def apply(localName: String): WomIdentifier = WomIdentifier(LocalName(localName), FullyQualifiedName(localName))
  def apply(localName: String, fullyQualifiedName: String): WomIdentifier =
    WomIdentifier(LocalName(localName), FullyQualifiedName(fullyQualifiedName))
}

case class WomIdentifier(localName: LocalName, fullyQualifiedName: FullyQualifiedName) {
  def combine(other: LocalName): WomIdentifier = combine(other.value)
  def combine(other: String): WomIdentifier =
    WomIdentifier(localName.combineToLocalName(other), fullyQualifiedName.combine(other))
  def workflowLocalName: String = fullyQualifiedName.value.split("\\.") match {
    case fqn if fqn.length > 1 => fqn.tail.mkString(".")
    case lqn => lqn.head
  }
}
