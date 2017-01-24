package cromwell.backend.sfs

import java.nio.file.Path

import cromwell.backend.standard.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}

class SharedFileSystemCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams)
  extends StandardCacheHitCopyingActor(standardParams) with SharedFileSystemJobCachingActorHelper {
  override protected def duplicate(source: Path, destination: Path): Unit = {
    // -Ywarn-value-discard
    sharedFileSystem.cacheCopy(source, destination)
    ()
  }
}
