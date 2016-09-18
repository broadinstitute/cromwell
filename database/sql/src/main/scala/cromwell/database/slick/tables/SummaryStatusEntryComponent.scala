package cromwell.database.slick.tables

import cromwell.database.sql.tables.SummaryStatus

trait SummaryStatusComponent {
  this: DriverComponent =>

  import driver.api._

  class SummaryStatuses(tag: Tag) extends Table[SummaryStatus](tag, "SUMMARY_STATUS") {
    def summaryStatusId = column[Int]("SUMMARY_STATUS_ID", O.PrimaryKey, O.AutoInc)

    def summaryTableName = column[String]("SUMMARY_TABLE_NAME")

    def summarizedTableName = column[String]("SUMMARIZED_TABLE_NAME")

    def maximumId = column[Long]("MAXIMUM_ID")

    override def * = (summaryTableName, summarizedTableName, maximumId, summaryStatusId.?) <>
      (SummaryStatus.tupled, SummaryStatus.unapply)

    def summarizedTableNameIndex = index("SUMMARY_STATUS_SUMMARY_TABLE_NAME_SUMMARIZED_TABLE_NAME_INDEX",
      (summaryTableName, summarizedTableName), unique = true)
  }

  protected val summaryStatuses = TableQuery[SummaryStatuses]

  val summaryStatusIdsAutoInc = summaryStatuses returning summaryStatuses.map(_.summaryStatusId)

  def summaryStatusMaximumIdBySummaryTableNameSummarizedTableName = Compiled(
    (summaryTableName: Rep[String], summarizedTableName: Rep[String]) => for {
      summaryStatus <- summaryStatuses
      if summaryStatus.summaryTableName === summaryTableName
      if summaryStatus.summarizedTableName === summarizedTableName
    } yield summaryStatus.maximumId
  )
}
