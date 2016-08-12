package cromwell.engine.workflow.lifecycle.execution

package object callcaching {

  // TODO: Find somewhere better for all these?
  case class HashKey(key: String)
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

  private[callcaching] case class CacheResultMatchesForHashes(hashResults: Iterable[HashResult], cacheResultIds: Set[Int])

  sealed trait CallCachingMode {
    def activity: Option[CallCachingActivity]

    def readFromCache = false
    def writeToCache = false
    def hashDockerNames = false
    def lookupDockerHashes: Boolean = false
    def hashFilePaths: Boolean = false
    def hashFileContents: Boolean = false

    /**
      * Return an equivalent of this call caching mode with READ disabled.
      * Especially useful to not re-check the cache again before the second attempt of a job.
      */
    def withoutRead: CallCachingMode
  }

  object CallCachingMode {
    def apply(read: Boolean, write: Boolean, hashDockerNames: Boolean, lookupDockerHashes: Boolean, hashFilePaths: Boolean, hashFileContents: Boolean): CallCachingMode = (read, write) match {
      case (false, false) => CallCachingOff
      case (false, true) => WriteCache(hashDockerNames, lookupDockerHashes, hashFilePaths, hashFileContents)
      case (true, false) => ReadCache(hashDockerNames, lookupDockerHashes, hashFilePaths, hashFileContents)
      case (true, true) => ReadAndWriteCache(hashDockerNames, lookupDockerHashes, hashFilePaths, hashFileContents)
    }
  }

  case object CallCachingOff extends CallCachingMode {
    override def activity = None
    override def withoutRead = this
  }

  sealed trait CallCachingActivity extends CallCachingMode {
    override def activity = Option(this)
  }
  case class ReadCache(override val hashDockerNames: Boolean, override val lookupDockerHashes: Boolean, override val hashFilePaths: Boolean, override val hashFileContents: Boolean) extends CallCachingActivity {
    override def readFromCache = true
    override def withoutRead = CallCachingOff
  }
  case class WriteCache(override val hashDockerNames: Boolean, override val lookupDockerHashes: Boolean, override val hashFilePaths: Boolean, override val hashFileContents: Boolean) extends CallCachingActivity {
    override def writeToCache = true
    override def withoutRead = this
  }
  case class ReadAndWriteCache(override val hashDockerNames: Boolean, override val lookupDockerHashes: Boolean, override val hashFilePaths: Boolean, override val hashFileContents: Boolean) extends CallCachingActivity {
    override def readFromCache = true
    override def writeToCache = true
    override def withoutRead = WriteCache(hashDockerNames, lookupDockerHashes, hashFilePaths, hashFileContents)
  }
}
