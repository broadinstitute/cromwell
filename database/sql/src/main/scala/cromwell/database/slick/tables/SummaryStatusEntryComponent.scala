package cromwell.database.slick.tables

import cromwell.database.sql.tables.SummaryStatusEntry

trait SummaryStatusEntryComponent {
  this: DriverComponent =>

  import driver.api._

  class SummaryStatusEntries(tag: Tag) extends Table[SummaryStatusEntry](tag, "SUMMARY_STATUS_ENTRY") {
    def summaryStatusEntryId = column[Int]("SUMMARY_STATUS_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def summaryTableName = column[String]("SUMMARY_TABLE_NAME")

    def summarizedTableName = column[String]("SUMMARIZED_TABLE_NAME")

    def maximumId = column[Long]("MAXIMUM_ID")

    override def * = (summaryTableName, summarizedTableName, maximumId, summaryStatusEntryId.?) <>
      (SummaryStatusEntry.tupled, SummaryStatusEntry.unapply)

    def ucSummaryStatusEntryStnStn = index("UC_SUMMARY_STATUS_ENTRY_STN_STN",
      (summaryTableName, summarizedTableName), unique = true)
  }

  protected val summaryStatusEntries = TableQuery[SummaryStatusEntries]

  val summaryStatusEntryIdsAutoInc = summaryStatusEntries returning summaryStatusEntries.map(_.summaryStatusEntryId)

  val maximumIdForSummaryTableNameSummarizedTableName = Compiled(
    (summaryTableName: Rep[String], summarizedTableName: Rep[String]) => for {
      summaryStatusEntry <- summaryStatusEntries
      if summaryStatusEntry.summaryTableName === summaryTableName
      if summaryStatusEntry.summarizedTableName === summarizedTableName
    } yield summaryStatusEntry.maximumId
  )
}
