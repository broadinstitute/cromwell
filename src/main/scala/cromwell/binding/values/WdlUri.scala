package cromwell.binding.values

import java.net.URI

import cromwell.binding.types.WdlUriType

import scala.util.{Success, Try}

// FIXME: THis will need to be linked to WdlFile (via a WdlFileLike trait or somesuch?) - see comment in https://broadinstitute.atlassian.net/browse/DSDEEPB-656)

object WdlUri {
  def apply(value: String): WdlUri = WdlUri(new URI(value))
}

case class WdlUri(value: URI) extends WdlPrimitive {
  val wdlType = WdlUriType

  override def add(rhs: WdlValue): Try[WdlValue] = ???

  override def equals(rhs: WdlValue): Try[WdlBoolean] = ???

  override def asString = value.toString
}
