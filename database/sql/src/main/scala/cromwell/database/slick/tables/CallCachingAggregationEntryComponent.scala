package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingAggregationEntry


trait CallCachingAggregationEntryComponent {

  this: DriverComponent with CallCachingEntryComponent with CallCachingDetritusEntryComponent =>

  import driver.api._

  class CallCachingAggregationEntries(tag: Tag) extends Table[CallCachingAggregationEntry](tag, "CALL_CACHING_AGGREGATION_ENTRY") {
    def callCachingAggregationEntryId = column[Int]("CALL_CACHING_AGGREGATION_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def baseAggregation = column[String]("BASE_AGGREGATION", O.Length(255))

    def inputFilesAggregation = column[Option[String]]("INPUT_FILES_AGGREGATION", O.Length(255))

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
  
  val callCachingAggregationForCacheEntryId = Compiled(
    (callCachingEntryId: Rep[Int]) => for {
      callCachingAggregationEntry <- callCachingAggregationEntries 
      if callCachingAggregationEntry.callCachingEntryId === callCachingEntryId
    } yield callCachingAggregationEntry
  )

  val existsCallCachingEntriesForBaseAggregationHash = Compiled(
    (baseAggregation: Rep[String]) => (for {
      callCachingEntry <- callCachingEntries
      if callCachingEntry.allowResultReuse
      callCachingAggregationEntry <- callCachingAggregationEntries
      if callCachingEntry.callCachingEntryId === callCachingAggregationEntry.callCachingEntryId
      if callCachingAggregationEntry.baseAggregation === baseAggregation
    } yield ()).exists
  )

  val existsCallCachingEntriesForBaseAggregationHashWithCallCachePrefix = Compiled(
    (baseAggregation: Rep[String],
     prefix1: Rep[String], prefix1Length: Rep[Int],
     prefix2: Rep[String], prefix2Length: Rep[Int],
     prefix3: Rep[String], prefix3Length: Rep[Int]
    ) => (for {
      callCachingEntry <- callCachingEntries
      if callCachingEntry.allowResultReuse
      callCachingAggregationEntry <- callCachingAggregationEntries
      if callCachingEntry.callCachingEntryId === callCachingAggregationEntry.callCachingEntryId
      if callCachingAggregationEntry.baseAggregation === baseAggregation
      detritus <- callCachingDetritusEntries
      if detritus.callCachingEntryId === callCachingEntry.callCachingEntryId
      detritusPath = detritus.detritusValue.map(clobToString)
      if (detritusPath.substring(0, prefix1Length) === prefix1) ||
        (detritusPath.substring(0, prefix2Length) === prefix2) ||
        (detritusPath.substring(0, prefix3Length) === prefix3)} yield ()).exists
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

  def callCachingEntriesForAggregatedHashesWithPrefixes(baseAggregation: Rep[String], inputFilesAggregation: Rep[Option[String]],
                                                        prefix1: Rep[String], prefix1Length: Rep[Int],
                                                        prefix2: Rep[String], prefix2Length: Rep[Int],
                                                        prefix3: Rep[String], prefix3Length: Rep[Int],
                                                        number: Int) = {
    (for {
      callCachingEntry <- callCachingEntries
      if callCachingEntry.allowResultReuse
      callCachingAggregationEntry <- callCachingAggregationEntries
      if callCachingEntry.callCachingEntryId === callCachingAggregationEntry.callCachingEntryId
      if callCachingAggregationEntry.baseAggregation === baseAggregation
      if (callCachingAggregationEntry.inputFilesAggregation.isEmpty && inputFilesAggregation.isEmpty) ||
        (callCachingAggregationEntry.inputFilesAggregation === inputFilesAggregation)
      detritus <- callCachingDetritusEntries
      // Pick only one detritus file since this is not an existence check and we don't want to return one row
      // for each of the (currently 6) types of standard detritus.
      if detritus.detritusKey === "returnCode"
      if detritus.callCachingEntryId === callCachingEntry.callCachingEntryId
      detritusPath = detritus.detritusValue.map(clobToString)
      if (detritusPath.substring(0, prefix1Length) === prefix1) ||
        (detritusPath.substring(0, prefix2Length) === prefix2) ||
        (detritusPath.substring(0, prefix3Length) === prefix3)
    } yield callCachingAggregationEntry.callCachingEntryId).drop(number - 1).take(1)
  }
}
