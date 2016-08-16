package cromwell.engine.workflow.lifecycle.execution

import cromwell.database.sql.MetaInfoId

package object callcaching {

  // TODO: Find somewhere better for all these?
  case class HashKey(key: String)
  object RuntimeAttributeHashKeyPlaceholder extends HashKey("PLACEHOLDER: Runtime Attributes")

  case class HashValue(value: String)
  case class HashResult(hashKey: HashKey, hashValue: HashValue)

  object HashValue {
    implicit class StringMd5er(unhashedString: String) {
      def md5HashValue: HashValue = {
        val hashBytes = java.security.MessageDigest.getInstance("MD5").digest(unhashedString.getBytes)
        HashValue(javax.xml.bind.DatatypeConverter.printHexBinary(hashBytes))
      }
    }
  }

  private[callcaching] case class CacheResultMatchesForHashes(hashResults: Set[HashResult], cacheResultIds: Set[MetaInfoId])
  private[callcaching] case class CacheResultLookupFailure(reason: Throwable)

  sealed trait CallCachingMode {
    /**
      * Return an equivalent of this call caching mode with READ disabled.
      */
    val withoutRead: CallCachingMode

    val readFromCache = false
    val writeToCache = false
  }

  case object CallCachingOff extends CallCachingMode {
    override val withoutRead = this
  }

  case class CallCachingActivity (readWriteMode: ReadWriteMode,
                                  dockerHashingType: DockerHashingType,
                                  fileHashingType: FileHashingType) extends CallCachingMode
  {
    override val readFromCache = readWriteMode.r
    override val writeToCache = readWriteMode.w
    override lazy val withoutRead: CallCachingMode = if (!writeToCache) CallCachingOff else this.copy(readWriteMode = WriteCache)
    override val toString = readWriteMode.toString
  }

  sealed trait ReadWriteMode {
    val r: Boolean = true
    val w: Boolean = true
  }
  case object ReadCache extends ReadWriteMode { override val w = false }
  case object WriteCache extends ReadWriteMode { override val r = false }
  case object ReadAndWriteCache extends ReadWriteMode

  sealed trait DockerHashingType
  case object HashDockerName extends DockerHashingType
  case object HashDockerNameAndLookupDockerHash extends DockerHashingType

  sealed trait FileHashingType
  case object HashFilePath extends FileHashingType
  case object HashFileContents extends FileHashingType
}
