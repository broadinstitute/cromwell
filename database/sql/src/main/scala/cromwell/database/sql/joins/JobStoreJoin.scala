package cromwell.database.sql.joins

import cromwell.database.sql.tables.{JobStoreEntry, JobStoreSimpletonEntry}

case class JobStoreJoin
(
  jobStoreEntry: JobStoreEntry,
  jobStoreSimpletonEntries: Seq[JobStoreSimpletonEntry]
)
