package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingHashEntry

trait CallCachingHashEntryComponent {

  this: DriverComponent with CallCachingEntryComponent =>

  import driver.api._

  class CallCachingHashEntries(tag: Tag) extends Table[CallCachingHashEntry](tag, "CALL_CACHING_HASH_ENTRY") {
    def callCachingHashEntryId = column[Long]("CALL_CACHING_HASH_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def hashKey = column[String]("HASH_KEY", O.Length(255))

    def hashValue = column[String]("HASH_VALUE", O.Length(255))

    def callCachingEntryId = column[Int]("CALL_CACHING_ENTRY_ID")

    override def * = (hashKey, hashValue, callCachingEntryId.?, callCachingHashEntryId.?) <>
      (CallCachingHashEntry.tupled, CallCachingHashEntry.unapply)

    def fkCallCachingHashEntryCallCachingEntryId = foreignKey("FK_CALL_CACHING_HASH_ENTRY_CALL_CACHING_ENTRY_ID",
      callCachingEntryId, callCachingEntries)(_.callCachingEntryId)

    def ucCallCachingHashEntryCceiHk =
      index("UC_CALL_CACHING_HASH_ENTRY_CCEI_HK", (callCachingEntryId, hashKey), unique = true)
  }

  val callCachingHashEntries = TableQuery[CallCachingHashEntries]

  val callCachingHashEntryIdsAutoInc = callCachingHashEntries returning
    callCachingHashEntries.map(_.callCachingHashEntryId)
  
  /**
    * Find all hashes for a CALL_CACHING_ENTRY_ID
    */
  val callCachingHashEntriesForCallCachingEntryId = Compiled(
    (callCachingEntryId: Rep[Int]) => for {
      callCachingHashEntry <- callCachingHashEntries
      if callCachingHashEntry.callCachingEntryId === callCachingEntryId
    } yield callCachingHashEntry
  )
}
