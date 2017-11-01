package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingAggregationEntry


trait CallCachingAggregationEntryComponent {

  this: DriverComponent with CallCachingEntryComponent =>

  import driver.api._

  class CallCachingAggregationEntries(tag: Tag) extends Table[CallCachingAggregationEntry](tag, "CALL_CACHING_AGGREGATION_ENTRY") {
    def callCachingAggregationEntryId = column[Int]("CALL_CACHING_AGGREGATION_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def baseAggregation = column[String]("BASE_AGGREGATION")

    def inputFilesAggregation = column[Option[String]]("INPUT_FILES_AGGREGATION")

    def callCachingEntryId = column[Int]("CALL_CACHING_ENTRY_ID")

    override def * = (baseAggregation, inputFilesAggregation, callCachingEntryId.?, callCachingAggregationEntryId.?) <>
      (CallCachingAggregationEntry.tupled, CallCachingAggregationEntry.unapply)

    def fkCallCachingAggregationEntryCallCachingEntryId = foreignKey("FK_CALL_CACHING_AGGREGATION_ENTRY_CALL_CACHING_ENTRY_ID",
      callCachingEntryId, callCachingEntries)(_.callCachingEntryId)

    def ixCallCachingAggregationEntryBaIfa =
      index("IX_CALL_CACHING_AGGREGATION_ENTRY_BA_IFA", (baseAggregation, inputFilesAggregation), unique = false)
  }

  val callCachingAggregationEntries = TableQuery[CallCachingAggregationEntries]

  val callCachingAggregationEntryIdsAutoInc = callCachingAggregationEntries returning
    callCachingAggregationEntries.map(_.callCachingAggregationEntryId)

  val existsCallCachingEntriesForBaseAggregationHash = Compiled(
    (baseAggregation: Rep[String]) => (for {
      callCachingEntry <- callCachingEntries
      if callCachingEntry.allowResultReuse
      callCachingAggregationEntry <- callCachingAggregationEntries
      if callCachingEntry.callCachingEntryId === callCachingAggregationEntry.callCachingEntryId
      if callCachingAggregationEntry.baseAggregation === baseAggregation
    } yield ()).exists
  )

  def callCachingEntriesForAggregatedHashes(baseAggregation: Rep[String], inputFilesAggregation: Rep[Option[String]], number: Int) = {
    (for {
      callCachingEntry <- callCachingEntries
      if callCachingEntry.allowResultReuse
      callCachingAggregationEntry <- callCachingAggregationEntries
      if callCachingEntry.callCachingEntryId === callCachingAggregationEntry.callCachingEntryId
      if callCachingAggregationEntry.baseAggregation === baseAggregation
      if (callCachingAggregationEntry.inputFilesAggregation.isEmpty && inputFilesAggregation.isEmpty) ||
        (callCachingAggregationEntry.inputFilesAggregation === inputFilesAggregation)
    } yield callCachingAggregationEntry.callCachingEntryId).drop(number - 1).take(1)
  }
}
