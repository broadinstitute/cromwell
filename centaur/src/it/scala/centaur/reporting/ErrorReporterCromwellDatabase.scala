package centaur.reporting

import cats.effect.IO
import centaur.CromwellDatabase
import cromwell.database.sql.tables.{JobKeyValueEntry, MetadataEntry}

import scala.concurrent.duration._
import scala.concurrent.ExecutionContext

class ErrorReporterCromwellDatabase(cromwellDatabase: CromwellDatabase) {

  def jobKeyValueEntriesIo(workflowExecutionUuidOption: Option[String])
                          (implicit executionContext: ExecutionContext): IO[Seq[JobKeyValueEntry]] = {
    workflowExecutionUuidOption.map(jobKeyValueEntriesIo).getOrElse(IO.pure(Seq.empty))
  }

  def jobKeyValueEntriesIo(workflowExecutionUuid: String)
                          (implicit executionContext: ExecutionContext): IO[Seq[JobKeyValueEntry]] = {
    IO.fromFuture(IO(cromwellDatabase.engineDatabase.queryJobKeyValueEntries(workflowExecutionUuid)))
  }

  def metadataEntriesIo(workflowExecutionUuidOption: Option[String])
                       (implicit executionContext: ExecutionContext): IO[Seq[MetadataEntry]] = {
    workflowExecutionUuidOption.map(metadataEntriesIo).getOrElse(IO.pure(Seq.empty))
  }

  def metadataEntriesIo(workflowExecutionUuid: String)
                       (implicit executionContext: ExecutionContext): IO[Seq[MetadataEntry]] = {
    // 30 seconds is less than production (60s as of 2018-08) but hopefully high enough to work on a CI machine with contended resources
    IO.fromFuture(IO(cromwellDatabase.metadataDatabase.queryMetadataEntries(workflowExecutionUuid, 30.seconds)))
  }

}
