package cromwell.database.slick.tables

import javax.sql.rowset.serial.SerialClob

import cromwell.database.sql.tables.CallCachingSimpletonEntry

trait CallCachingSimpletonEntryComponent {

  this: DriverComponent with CallCachingEntryComponent =>

  import driver.api._

  class CallCachingSimpletonEntries(tag: Tag)
    extends Table[CallCachingSimpletonEntry](tag, "CALL_CACHING_SIMPLETON_ENTRY") {
    def callCachingSimpletonEntryId = column[Int]("CALL_CACHING_SIMPLETON_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def simpletonKey = column[String]("SIMPLETON_KEY", O.Length(255))

    def simpletonValue = column[Option[SerialClob]]("SIMPLETON_VALUE")

    def wdlType = column[String]("WDL_TYPE", O.Length(255))

    def callCachingEntryId = column[Int]("CALL_CACHING_ENTRY_ID")

    override def * = (simpletonKey, simpletonValue, wdlType, callCachingEntryId.?, callCachingSimpletonEntryId.?) <>
      (CallCachingSimpletonEntry.tupled, CallCachingSimpletonEntry.unapply)

    def fkCallCachingSimpletonEntryCallCachingEntryId = foreignKey(
      "FK_CALL_CACHING_SIMPLETON_ENTRY_CALL_CACHING_ENTRY_ID", callCachingEntryId, callCachingEntries)(_.callCachingEntryId)

    def ucCallCachingSimpletonEntryCceiSk =
      index("UC_CALL_CACHING_SIMPLETON_ENTRY_CCEI_SK", (callCachingEntryId, simpletonKey), unique = true)
  }

  val callCachingSimpletonEntries = TableQuery[CallCachingSimpletonEntries]

  val callCachingSimpletonEntryIdsAutoInc = callCachingSimpletonEntries returning
    callCachingSimpletonEntries.map(_.callCachingSimpletonEntryId)

  /**
    * Find all result simpletons which match a given CALL_CACHING_ENTRY_ID
    */
  val callCachingSimpletonEntriesForCallCachingEntryId = Compiled(
    (callCachingEntryId: Rep[Int]) => for {
      callCachingSimpletonEntry <- callCachingSimpletonEntries
      if callCachingSimpletonEntry.callCachingEntryId === callCachingEntryId
    } yield callCachingSimpletonEntry
  )
}
