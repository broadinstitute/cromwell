package cromwell.binding.values

import cromwell.engine._

trait WdlPrimitive extends WdlValue {

  /**
   * No need to hash primitives.
   */
  override def getHash(implicit hasher: FileHasher) = valueString
}
