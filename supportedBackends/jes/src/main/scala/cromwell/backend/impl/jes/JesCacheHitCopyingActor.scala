package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.services.io.AsyncIo

import scala.concurrent.Future

import cromwell.backend.standard.{StandardCacheHitCopyingActor, StandardCacheHitCopyingActorParams}

case class JesCacheHitCopyingActor(override val standardParams: StandardCacheHitCopyingActorParams)
  extends StandardCacheHitCopyingActor with JesJobCachingActorHelper with AsyncIo {
  override protected def duplicate(source: Path, destination: Path): Future[Unit] = copy(source, destination)
}
