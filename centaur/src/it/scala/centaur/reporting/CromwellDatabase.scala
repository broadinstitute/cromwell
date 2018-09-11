package centaur.reporting

import cats.effect.IO
import cromwell.database.sql.tables.{JobKeyValueEntry, MetadataEntry}
import cromwell.database.sql.{EngineSqlDatabase, MetadataSqlDatabase}

import scala.concurrent.ExecutionContext

/**
  * Wraps connections to a cromwell database. The database connections are not initialized until first use.
  */
class CromwellDatabase(engineDatabaseThunk: => EngineSqlDatabase, metadataDatabaseThunk: => MetadataSqlDatabase) {

  private lazy val engineDatabase: EngineSqlDatabase = engineDatabaseThunk
  private lazy val metadataDatabase: MetadataSqlDatabase = metadataDatabaseThunk

  def jobKeyValueEntriesIo(workflowExecutionUuidOption: Option[String])
                          (implicit executionContext: ExecutionContext): IO[Seq[JobKeyValueEntry]] = {
    workflowExecutionUuidOption.map(jobKeyValueEntriesIo).getOrElse(IO.pure(Seq.empty))
  }

  def jobKeyValueEntriesIo(workflowExecutionUuid: String)
                          (implicit executionContext: ExecutionContext): IO[Seq[JobKeyValueEntry]] = {
    IO.fromFuture(IO(engineDatabase.queryJobKeyValueEntries(workflowExecutionUuid)))
  }

  def metadataEntriesIo(workflowExecutionUuidOption: Option[String])
                       (implicit executionContext: ExecutionContext): IO[Seq[MetadataEntry]] = {
    workflowExecutionUuidOption.map(metadataEntriesIo).getOrElse(IO.pure(Seq.empty))
  }

  def metadataEntriesIo(workflowExecutionUuid: String)
                       (implicit executionContext: ExecutionContext): IO[Seq[MetadataEntry]] = {
    IO.fromFuture(IO(metadataDatabase.queryMetadataEntries(workflowExecutionUuid)))
  }

}
