package cromwell.binding.values

import cromwell.engine._

trait WdlPrimitive extends WdlValue {

  /**
   * No need to hash primitives but add the className to the value to differentiate WdlTypes.
   */
  override def getHash(implicit hasher: FileHasher) = getClass.getCanonicalName+valueString
}
