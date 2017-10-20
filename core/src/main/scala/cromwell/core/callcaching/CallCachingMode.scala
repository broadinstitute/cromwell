package cromwell.core.callcaching

sealed trait CallCachingMode {
  /**
    * Return an equivalent of this call caching mode with READ disabled.
    */
  val withoutRead: CallCachingMode

  val withoutWrite: CallCachingMode

  val readFromCache = false
  val writeToCache = false
}

case object CallCachingOff extends CallCachingMode {
  override val readFromCache = false
  override val writeToCache = false
  override val withoutRead = this
  override val withoutWrite = this
}

case class CallCachingActivity(readWriteMode: ReadWriteMode, options: CallCachingOptions = CallCachingOptions(invalidateBadCacheResults = true)) extends CallCachingMode {
  override val readFromCache = readWriteMode.r
  override val writeToCache = readWriteMode.w
  override lazy val withoutRead: CallCachingMode = if (!writeToCache) CallCachingOff else this.copy(readWriteMode = WriteCache)
  override lazy val withoutWrite: CallCachingMode = if (!readFromCache) CallCachingOff else this.copy(readWriteMode = ReadCache)
  override val toString = readWriteMode.toString
}

sealed trait ReadWriteMode {
  val r: Boolean = true
  val w: Boolean = true
}
case object ReadCache extends ReadWriteMode { override val w = false }
case object WriteCache extends ReadWriteMode { override val r = false }
case object ReadAndWriteCache extends ReadWriteMode

final case class CallCachingOptions(invalidateBadCacheResults: Boolean = true)
