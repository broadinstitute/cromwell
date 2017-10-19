package cromwell

import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.tables.{CallCachingSimpletonEntry, JobStoreSimpletonEntry}
import wom.types.{WdlBooleanType, WdlFloatType, WdlIntegerType, WdlStringType}
import wom.values.{WdlPrimitive, WdlSingleFile, WdlValue}

import scala.util.Try

/** Converts database simpletons to `WdlValueSimpleton`s, explicitly instantiating objects of "File" type
  * to `WdlSingleFile` instances.
  */
object Simpletons {
  def toSimpleton(entry: CallCachingSimpletonEntry): WdlValueSimpleton = {
    toSimpleton(entry.wdlType, entry.simpletonKey, entry.simpletonValue.toRawString)
  }

  def toSimpleton(entry: JobStoreSimpletonEntry): WdlValueSimpleton = {
    toSimpleton(entry.wdlType, entry.simpletonKey, entry.simpletonValue.toRawString)
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
