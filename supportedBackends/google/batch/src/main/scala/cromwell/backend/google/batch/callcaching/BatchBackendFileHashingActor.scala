package cromwell.backend.google.batch.callcaching

import cromwell.backend.standard.callcaching.{
  AsyncFileHashingStrategy,
  StandardFileHashingActor,
  StandardFileHashingActorParams
}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder

class BatchBackendFileHashingActor(standardParams: StandardFileHashingActorParams)
    extends StandardFileHashingActor(standardParams) {
  override val ioCommandBuilder = GcsBatchCommandBuilder

  override val defaultHashingStrategies: Map[String, AsyncFileHashingStrategy] = Map(
    ("gcs", AsyncFileHashingStrategy.Crc32c),
    ("drs", AsyncFileHashingStrategy.Crc32c)
  )
}
