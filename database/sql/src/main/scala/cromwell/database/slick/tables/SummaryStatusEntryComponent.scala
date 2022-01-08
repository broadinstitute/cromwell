package cromwell.database.slick.tables

import cromwell.database.sql.tables.SummaryStatusEntry

trait SummaryStatusEntryComponent {
  this: DriverComponent =>

  import driver.api._

  class SummaryStatusEntries(tag: Tag) extends Table[SummaryStatusEntry](tag, "SUMMARY_STATUS_ENTRY") {
    def summaryStatusEntryId = column[Int]("SUMMARY_STATUS_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def summaryName = column[String]("SUMMARY_NAME", O.Length(255))

    def summaryPosition = column[Long]("SUMMARY_POSITION")

    override def * = (summaryName, summaryPosition, summaryStatusEntryId.?) <>
      (SummaryStatusEntry.tupled, SummaryStatusEntry.unapply)

    def ucSummaryStatusEntrySn = index("UC_SUMMARY_STATUS_ENTRY_SN", summaryName, unique = true)
  }

  protected val summaryStatusEntries = TableQuery[SummaryStatusEntries]

  val summaryStatusEntryIdsAutoInc = summaryStatusEntries returning summaryStatusEntries.map(_.summaryStatusEntryId)

  val summaryPositionForSummaryName = Compiled(
    (summaryName: Rep[String]) => for {
      summaryStatusEntry <- summaryStatusEntries
      if summaryStatusEntry.summaryName === summaryName
    } yield summaryStatusEntry.summaryPosition
  )
}
