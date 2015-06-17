package cromwell.binding.values

import java.nio.file.{Paths, Path}
import cromwell.binding.types.WdlFileType

object WdlFile {
  def apply(value: String): WdlFile = WdlFile(Paths.get(value))
}

import scala.util.{Success, Try}

case class WdlFile(value: Path) extends WdlPrimitive {
  val wdlType = WdlFileType
  override def add(rhs: WdlValue): Try[WdlValue] = {
    rhs match {
      // The goofiness with the value.toStrings are because we want concatenation and not Path.resolve() logic
      case r: WdlString => Success(WdlFile(Paths.get(value.toString, r.value)))
      case r: WdlFile => Success(WdlFile(Paths.get(value.toString, r.value.toString)))
      case _ => invalid(s"$value + $rhs")
    }
  }
  override def equals(rhs: WdlValue): Try[WdlBoolean] = {
    rhs match {
      case r: WdlFile => Success(WdlBoolean(value.equals(r.value)))
      case r: WdlString => Success(WdlBoolean(value.toString.equals(r.value)))
      case _ => invalid(s"$value == $rhs")
    }
  }
  override def toWdlString = value.toString
}
