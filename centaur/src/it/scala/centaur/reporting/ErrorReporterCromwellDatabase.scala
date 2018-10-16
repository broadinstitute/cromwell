package centaur.reporting

import cats.effect.IO
import centaur.CromwellDatabase
import cromwell.database.sql.tables.{JobKeyValueEntry, MetadataEntry}

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
    locally(executionContext)
//    IO.fromFuture(IO(cromwellDatabase.metadataDatabase.queryMetadataEntries(workflowExecutionUuid)))
    IO.pure(List.empty)
  }

}
