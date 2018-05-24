package wdl.draft3.transforms.wdlom2wdl

import simulacrum.typeclass
import scala.language.implicitConversions

@typeclass
trait WdlWriter[A] {
  def toWdlV1(a: A): String
}

object WdlWriter {
  // Stolen from WomGraph.scala
  final def indent(s: String) = s.lines.map(x => s"  $x").mkString(System.lineSeparator)
  final def combine(ss: Iterable[String]) = ss.mkString(start="", sep=System.lineSeparator, end=System.lineSeparator)
  final def indentAndCombine(ss: Iterable[String]) = combine(ss.map(indent))
}
