package cromwell.core.callcaching

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
                                dockerHashingType: DockerHashingType = HashDockerName) extends CallCachingMode
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
