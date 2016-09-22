package cromwell.database.slick.tables

import cromwell.database.sql.tables.CallCachingDetritusEntry

trait CallCachingDetritusEntryComponent {

  this: DriverComponent with CallCachingEntryComponent =>

  import driver.api._

  class CallCachingDetritusEntries(tag: Tag)
    extends Table[CallCachingDetritusEntry](tag, "CALL_CACHING_DETRITUS_ENTRY") {
    def callCachingDetritusEntryId = column[Int]("CALL_CACHING_DETRITUS_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def detritusKey = column[String]("DETRITUS_KEY")

    def detritusValue = column[String]("DETRITUS_VALUE")

    def callCachingEntryId = column[Int]("CALL_CACHING_ENTRY_ID")

    override def * = (detritusKey, detritusValue, callCachingEntryId.?, callCachingDetritusEntryId.?) <>
      (CallCachingDetritusEntry.tupled, CallCachingDetritusEntry.unapply)

    def fkCallCachingDetritusEntryCallCachingEntryId = foreignKey(
      "FK_CALL_CACHING_DETRITUS_ENTRY_CALL_CACHING_ENTRY_ID",
      callCachingEntryId, callCachingEntries)(_.callCachingEntryId)

    def ucCallCachingDetritusEntryCceiDk =
      index("UC_CALL_CACHING_DETRITUS_ENTRY_CCEI_DK", (callCachingEntryId, detritusKey), unique = true)
  }

  protected val callCachingDetritusEntries = TableQuery[CallCachingDetritusEntries]

  val callCachingDetritusEntryIdsAutoInc =
    callCachingDetritusEntries returning callCachingDetritusEntries.map(_.callCachingDetritusEntryId)

  val callCachingDetritusEntriesForCallCachingEntryId = Compiled(
    (callCachingEntryId: Rep[Int]) => for {
      callCachingDetritusEntry <- callCachingDetritusEntries
      if callCachingDetritusEntry.callCachingEntryId === callCachingEntryId
    } yield callCachingDetritusEntry
  )
}
