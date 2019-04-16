package cromwell.backend.impl.bcs.callcaching


import cromwell.backend.standard.callcaching.{StandardFileHashingActor, StandardFileHashingActorParams}
import cromwell.filesystems.oss.batch.OssBatchCommandBuilder

class BcsBackendFileHashingActor(standardParams: StandardFileHashingActorParams) extends StandardFileHashingActor(standardParams) {
  override val ioCommandBuilder = OssBatchCommandBuilder
}