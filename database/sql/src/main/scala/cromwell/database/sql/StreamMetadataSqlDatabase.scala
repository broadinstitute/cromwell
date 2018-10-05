package cromwell.database.sql

import cats.data.NonEmptyList
import cromwell.database.sql.joins.MetadataJobQueryValue
import cromwell.database.sql.tables.MetadataEntry
import slick.basic.DatabasePublisher

trait StreamMetadataSqlDatabase extends SqlDatabase {

  def streamQueryMetadataEntries(workflowExecutionUuid: String): DatabasePublisher[MetadataEntry]

  def streamQueryMetadataEntries(workflowExecutionUuid: String,
                           metadataKey: String): DatabasePublisher[MetadataEntry]

  def streamQueryMetadataEntries(workflowExecutionUuid: String,
                           callFullyQualifiedName: String,
                           jobIndex: Option[Int],
                           jobAttempt: Option[Int]): DatabasePublisher[MetadataEntry]

  def streamQueryMetadataEntries(workflowUuid: String,
                           metadataKey: String,
                           callFullyQualifiedName: String,
                           jobIndex: Option[Int],
                           jobAttempt: Option[Int]): DatabasePublisher[MetadataEntry]

  def streamQueryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String,
                                           metadataKeys: NonEmptyList[String],
                                           metadataJobQueryValue: MetadataJobQueryValue): DatabasePublisher[MetadataEntry]

  def streamQueryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String,
                                            metadataKeys: NonEmptyList[String],
                                            metadataJobQueryValue: MetadataJobQueryValue): DatabasePublisher[MetadataEntry]
}
