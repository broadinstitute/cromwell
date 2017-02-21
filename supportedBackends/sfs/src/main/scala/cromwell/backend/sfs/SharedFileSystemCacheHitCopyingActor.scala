package cromwell.backend.sfs

import cromwell.backend.standard.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.core.path.Path

class SharedFileSystemCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams)
  extends StandardCacheHitCopyingActor(standardParams) with SharedFileSystemJobCachingActorHelper {
  override protected def duplicate(source: Path, destination: Path): Unit = {
    sharedFileSystem.cacheCopy(source, destination).get
  }
}
