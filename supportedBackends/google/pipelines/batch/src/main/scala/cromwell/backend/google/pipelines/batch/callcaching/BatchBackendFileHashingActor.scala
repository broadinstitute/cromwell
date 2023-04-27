package cromwell.backend.google.pipelines.batch.callcaching

import cromwell.backend.standard.callcaching.{StandardFileHashingActor, StandardFileHashingActorParams}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder

class BatchBackendFileHashingActor(standardParams: StandardFileHashingActorParams) extends StandardFileHashingActor(standardParams) {
  override val ioCommandBuilder = GcsBatchCommandBuilder
}
