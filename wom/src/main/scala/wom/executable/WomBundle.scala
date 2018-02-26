package wom.executable

import wom.callable.Callable
import wom.types.WomType

/**
  * Represents a set of static WOM items that might be imported from a file.
  */
final case class WomBundle(callables: Set[Callable],
                           typeAliases: Map[String, WomType])
