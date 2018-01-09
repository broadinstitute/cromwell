package cromwell.backend.impl.jes.callcaching

import cromwell.backend.standard.callcaching.{StandardFileHashingActor, StandardFileHashingActorParams}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder

class JesBackendFileHashingActor(standardParams: StandardFileHashingActorParams) extends StandardFileHashingActor(standardParams) {
  override val ioCommandBuilder = GcsBatchCommandBuilder
}
