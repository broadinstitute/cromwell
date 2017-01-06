package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.standard.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}
import cromwell.core.path.PathCopier

case class JesCacheHitCopyingActor(override val standardParams: StandardCacheHitCopyingActorParams)
  extends StandardCacheHitCopyingActor with JesJobCachingActorHelper {
  override protected def duplicate(source: Path, destination: Path): Unit = PathCopier.copy(source, destination).get
}
