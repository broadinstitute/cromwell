package cromwell

import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.database.sql.tables.{CallCachingSimpletonEntry, JobStoreSimpletonEntry}
import wdl4s.types.{WdlBooleanType, WdlFloatType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlPrimitive, WdlSingleFile, WdlValue}

import scala.util.Try

/** Converts database simpletons to `WdlValueSimpleton`s, explicitly instantiating objects of "File" type
  * to `WdlSingleFile` instances.
  */
object Simpletons {
  def toSimpleton(entry: CallCachingSimpletonEntry): WdlValueSimpleton = {
    toSimpleton(entry.wdlType, entry.simpletonKey, entry.simpletonValue)
  }

  def toSimpleton(entry: JobStoreSimpletonEntry): WdlValueSimpleton = {
    toSimpleton(entry.wdlType, entry.simpletonKey, entry.simpletonValue)
  }

  private def toSimpleton(wdlType: String, simpletonKey: String, simpletonValue: String): WdlValueSimpleton = {
    val wdlValue: String => Try[WdlValue] = wdlType match {
      case "String" => WdlStringType.coerceRawValue
      case "Int" => WdlIntegerType.coerceRawValue
      case "Float" => WdlFloatType.coerceRawValue
      case "Boolean" => WdlBooleanType.coerceRawValue
      case "File" => s => Try(WdlSingleFile(s))
      case _ => throw new RuntimeException(s"unrecognized WDL type: $wdlType")
    }
    WdlValueSimpleton(simpletonKey, wdlValue(simpletonValue).get.asInstanceOf[WdlPrimitive])
  }
}
