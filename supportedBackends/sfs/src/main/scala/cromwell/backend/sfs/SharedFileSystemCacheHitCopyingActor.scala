package cromwell.backend.sfs

import cromwell.backend.standard.callcaching.StandardCacheHitCopyingActor.PathPair
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.filesystems.gcs.GcsBatchCommandBuilder
import lenthall.util.TryUtil

import scala.util.{Failure, Try}

class SharedFileSystemCacheHitCopyingActor(standardParams: StandardCacheHitCopyingActorParams)
  extends StandardCacheHitCopyingActor(standardParams) with SharedFileSystemJobCachingActorHelper with GcsBatchCommandBuilder {
  override protected def duplicate(copyPairs: Set[PathPair]): Option[Try[Unit]] = Option {
    val copies = copyPairs map {
      case (source, destination) => 
        sharedFileSystem.cacheCopy(source, destination)
    }
    
    TryUtil.sequence(copies.toList) map { _ => () } recoverWith {
      case failure =>
        // If one or more of the copies failed, we want to delete all the files that were successfully copied
        // before that. Especially if they've been symlinked, leaving them could lead to rewriting the original
        // files when the job gets re-run
        // TODO: this could be done more generally in the StandardCacheHitCopyingActor
        copyPairs foreach {
          case (_, dst) => dst.delete(swallowIOExceptions = true)
        }
        Failure(failure)
    }
  }
}
