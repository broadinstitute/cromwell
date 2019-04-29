// tl;dr postgresql + jdbc + slick = :(

package cromwell.database.sql

import java.sql.{Blob, Clob}
import javax.sql.rowset.serial.{SerialBlob, SerialClob}

import cromwell.database.sql.tables._
import cromwell.database.sql.joins.{CallCachingJoin, JobStoreJoin}

/* XXX Postgres large object workaround
 * These implicits cope with an irritating behavior of the Postgres JDBC
 * driver, which appears to lazy-load Blob and Clob fields.  By immediately
 * converting to javax.sql types within a transaction, we prevent future
 * database calls.  This should not be used unless Cromwell is configured for
 * Postgresql; MySql stores large objects directly in the tables.
 *
 * Note that only Blob is actually using the 'lo' extension in the tables,
 * Clob is stored as text, but it still triggers the transaction error.
 */
object SqlTableConverters {

  private def toSerialClob(clob: Clob): SerialClob = new SerialClob(clob)
  private def toSerialBlob(blob: Blob): SerialBlob = new SerialBlob(blob)

  implicit class WorkflowStoreEntryTransactional(val entry: WorkflowStoreEntry)
      extends AnyVal {
    def withLargeObjects = entry.copy(
      workflowDefinition = entry.workflowDefinition.map(toSerialClob),
      workflowInputs = entry.workflowInputs.map(toSerialClob),
      workflowOptions = entry.workflowOptions.map(toSerialClob),
      importsZip = entry.importsZip.map(toSerialBlob),
      customLabels = toSerialClob(entry.customLabels))
  }

  implicit class MetadataEntryTransactional(val entry: MetadataEntry)
      extends AnyVal {
    def withLargeObjects = entry.copy(
      metadataValue = entry.metadataValue.map(toSerialClob))
  }

  implicit class JobStoreEntryTransactional(val entry: JobStoreEntry)
      extends AnyVal {
    def withLargeObjects = entry.copy(
      exceptionMessage = entry.exceptionMessage.map(toSerialClob))
  }

  implicit class JobStoreSimpletonEntryTransactional(val entry: JobStoreSimpletonEntry)
      extends AnyVal {
    def withLargeObjects = entry.copy(
      simpletonValue = entry.simpletonValue.map(toSerialClob))
  }

  implicit class CallCachingSimpletonEntryTransactional(val entry: CallCachingSimpletonEntry)
      extends AnyVal {
    def withLargeObjects = entry.copy(
      simpletonValue = entry.simpletonValue.map(toSerialClob))
  }

  implicit class CallCachingDetritusEntryTransactional(val entry: CallCachingDetritusEntry)
      extends AnyVal {
    def withLargeObjects = entry.copy(
      detritusValue = entry.detritusValue.map(toSerialClob))
  }

  implicit class CallCachingJoinTransactional(val join: CallCachingJoin)
      extends AnyVal {
    def withLargeObjects = join.copy(
      callCachingSimpletonEntries = join.callCachingSimpletonEntries.map(
        _.withLargeObjects),
      callCachingDetritusEntries = join.callCachingDetritusEntries.map(
        _.withLargeObjects))
  }

  implicit class JobStoreJoinTransactional(val join: JobStoreJoin)
      extends AnyVal {
    def withLargeObjects = join.copy(
      jobStoreEntry = join.jobStoreEntry.withLargeObjects,
      jobStoreSimpletonEntries = join.jobStoreSimpletonEntries.map(_.withLargeObjects))
  }
}
