package cromwell.binding.values

import cromwell.binding.types.WdlFileType

/**
 * WDL value which identifies the location of a File-like object.
 */
trait WdlFileLike extends WdlPrimitive {
  val wdlType = WdlFileType
}
