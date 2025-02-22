package cromwell.backend.google.batch.callcaching

import cromwell.backend.standard.callcaching.{StandardFileHashingActor, StandardFileHashingActorParams}
import cromwell.core.callcaching.FileHashStrategy
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder

class BatchBackendFileHashingActor(standardParams: StandardFileHashingActorParams)
    extends StandardFileHashingActor(standardParams) {
  override val ioCommandBuilder = GcsBatchCommandBuilder(metricsCallback)

  override val defaultHashingStrategies: Map[String, FileHashStrategy] = Map(
    ("gcs", FileHashStrategy.Crc32c),
    ("drs", FileHashStrategy.Drs)
  )
}
