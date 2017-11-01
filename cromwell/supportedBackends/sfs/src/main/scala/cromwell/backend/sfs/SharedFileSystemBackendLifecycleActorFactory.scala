package cromwell.backend.sfs

import cromwell.backend.standard._
import cromwell.backend.standard.callcaching.StandardCacheHitCopyingActor

/**
  * A factory that can be extended for any shared file system implementation.
  *
  * See the SharedFileSystemAsyncJobExecutionActor for more info.
  */
trait SharedFileSystemBackendLifecycleActorFactory extends StandardLifecycleActorFactory {

  override def jobIdKey: String = SharedFileSystemAsyncJobExecutionActor.JobIdKey

  override lazy val cacheHitCopyingActorClassOption: Option[Class[_ <: StandardCacheHitCopyingActor]] = {
    Option(classOf[SharedFileSystemCacheHitCopyingActor])
  }
}
