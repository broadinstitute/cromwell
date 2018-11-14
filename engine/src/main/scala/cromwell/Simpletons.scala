package cromwell

import cromwell.core.simpleton.WomValueSimpleton
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.tables.{CallCachingSimpletonEntry, JobStoreSimpletonEntry}
import wom.types.{WomBooleanType, WomFloatType, WomIntegerType, WomLongType, WomStringType}
import wom.values.{WomPrimitive, WomSingleFile, WomUnlistedDirectory, WomValue}

import scala.util.Try

/** Converts database simpletons to `WomValueSimpleton`s, explicitly instantiating objects of "File" type
  * to `WdlSingleFile` instances.
  */
object Simpletons {
  def toSimpleton(entry: CallCachingSimpletonEntry): WomValueSimpleton = {
    toSimpleton(entry.wdlType, entry.simpletonKey, entry.simpletonValue.toRawString)
  }

  def toSimpleton(entry: JobStoreSimpletonEntry): WomValueSimpleton = {
    toSimpleton(entry.wdlType, entry.simpletonKey, entry.simpletonValue.toRawString)
  }

  private def toSimpleton(womType: String, simpletonKey: String, simpletonValue: String): WomValueSimpleton = {
    val womValue: String => Try[WomValue] = womType match {
      case "String" => WomStringType.coerceRawValue
      case "Int" => WomIntegerType.coerceRawValue
      case "Float" => WomFloatType.coerceRawValue
      case "Long" => WomLongType.coerceRawValue
      case "Boolean" => WomBooleanType.coerceRawValue
      case "File" => s => Try(WomSingleFile(s))
      case "Directory" => s => Try(WomUnlistedDirectory(s))
      case _ => throw new RuntimeException(s"unrecognized simpleton WOM type: $womType")
    }
    WomValueSimpleton(simpletonKey, womValue(simpletonValue).get.asInstanceOf[WomPrimitive])
  }
}
