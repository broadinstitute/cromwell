package cromwell

import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.database.sql.tables.DatabaseSimpleton
import wdl4s.types.{WdlBooleanType, WdlFloatType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlPrimitive, WdlSingleFile, WdlValue}

import scala.util.Try

object Simpletons {
  /** Converts `DatabaseSimpleton`s to `WdlValueSimpleton`s, explicitly instantiating objects of "File" type
    * to `WdlSingleFile` instances.
    */
  def toSimpleton(entry: DatabaseSimpleton): WdlValueSimpleton = {
    val wdlValue: String => Try[WdlValue] = entry.wdlType match {
      case "String" => WdlStringType.coerceRawValue
      case "Int" => WdlIntegerType.coerceRawValue
      case "Float" => WdlFloatType.coerceRawValue
      case "Boolean" => WdlBooleanType.coerceRawValue
      case "File" => s => Try(WdlSingleFile(s))
      case _ => throw new RuntimeException(s"$entry: unrecognized WDL type: ${entry.wdlType}")
    }
    WdlValueSimpleton(entry.simpletonKey, wdlValue(entry.simpletonValue).get.asInstanceOf[WdlPrimitive])
  }
}
