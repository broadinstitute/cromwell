package cromwell.engine.workflow.lifecycle.execution

package object callcaching {

  // TODO: Find somewhere better for all these?
  case class HashKey(key: String)
  case class HashValue(value: String)
  case class HashResult(hashKey: HashKey, hashValue: HashValue)

  private[callcaching] case class CacheResultMatchesForHashes(hashResults: Iterable[HashResult], cacheResultIds: Set[Int])

  sealed trait CallCachingMode { def activity: Option[CallCachingActivity]; def readFromCache = false; def writeToCache = false; def lookupDockerHashes: Boolean = false }
  object CallCachingMode {
    def apply(read: Boolean, write: Boolean, lookupDockerHashes: Boolean): CallCachingMode = (read, write) match {
      case (false, false) => CallCachingOff
      case (false, true) => WriteCache(lookupDockerHashes)
      case (true, false) => ReadCache(lookupDockerHashes)
      case (true, true) => ReadAndWriteCache(lookupDockerHashes)
    }
  }

  case object CallCachingOff extends CallCachingMode { override def activity = None; }

  sealed trait CallCachingActivity extends CallCachingMode { override def activity = Option(this); }
  case class ReadCache(override val lookupDockerHashes: Boolean) extends CallCachingActivity { override def readFromCache = true }
  case class WriteCache(override val lookupDockerHashes: Boolean) extends CallCachingActivity { override def writeToCache = true }
  case class ReadAndWriteCache(override val lookupDockerHashes: Boolean) extends CallCachingActivity { override def readFromCache = true; override def writeToCache = true }
}
