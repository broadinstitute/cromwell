package cromwell.backend.sfs

import java.nio.file.Path

import akka.actor.ActorRef
import cromwell.backend.callcaching.CacheHitDuplicating
import cromwell.backend.{BackendCacheHitCopyingActor, BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptor}
import cromwell.core.path.PathFactory

import scala.util.Try

class SharedFileSystemCacheHitCopyingActor(override val jobDescriptor: BackendJobDescriptor,
                                           override val configurationDescriptor: BackendConfigurationDescriptor,
                                           override val backendInitializationDataOption:
                                           Option[BackendInitializationData],
                                           override val serviceRegistryActor: ActorRef)
  extends SharedFileSystemJobCachingActorHelper with BackendCacheHitCopyingActor with CacheHitDuplicating {

  override lazy val destinationCallRootPath = jobPaths.callRoot

  override lazy val destinationJobDetritusPaths = jobPaths.detritusPaths

  override protected def getPath(file: String) = Try(PathFactory.buildPath(file, jobPaths.pathBuilders))

  override protected def duplicate(source: Path, destination: Path) = {
    // -Ywarn-value-discard
    sharedFileSystem.cacheCopy(source, destination)
    ()
  }
}
