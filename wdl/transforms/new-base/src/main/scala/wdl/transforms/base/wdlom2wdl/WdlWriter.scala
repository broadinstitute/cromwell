package wdl.transforms.base.wdlom2wdl

import simulacrum.typeclass
import scala.language.implicitConversions

@typeclass
trait WdlWriter[A] {
  def toWdlV1(a: A): String
}

object WdlWriter {
  // Stolen from WomGraph.scala
  def indent(s: String) = s.linesIterator.map(x => s"  $x").mkString(System.lineSeparator)
  def combine(ss: Iterable[String]) = ss.mkString(start="", sep=System.lineSeparator, end=System.lineSeparator)
  def indentAndCombine(ss: Iterable[String]) = combine(ss.map(indent))
}
