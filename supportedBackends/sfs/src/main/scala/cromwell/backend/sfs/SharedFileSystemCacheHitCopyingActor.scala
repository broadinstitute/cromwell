package cromwell.backend.sfs

import java.nio.file.Path

import cromwell.backend.standard.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}

import scala.concurrent.Future

class SharedFileSystemCacheHitCopyingActor(override val standardParams: StandardCacheHitCopyingActorParams)
  extends StandardCacheHitCopyingActor with SharedFileSystemJobCachingActorHelper {
  override protected def duplicate(source: Path, destination: Path): Future[Unit] = {
    // -Ywarn-value-discard
    sharedFileSystem.cacheCopy(source, destination)
    Future.successful(())
  }
}
