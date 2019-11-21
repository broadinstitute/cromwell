package cromwell.database.slick

import scala.concurrent.ExecutionContext

trait SummaryStatusSlickDatabase {
  this: MetadataSlickDatabase =>

  import dataAccess.driver.api._

  private[slick] def getSummaryStatusEntrySummaryPosition(summaryName: String): DBIO[Option[Long]] = {
    dataAccess.summaryPositionForSummaryName(summaryName).result.headOption
  }

  private[slick] def upsertSummaryStatusEntrySummaryPosition(summaryName: String,
                                                             oldSummaryPosition: Long,
                                                             newSummaryPosition: Long)
                                                            (implicit ec: ExecutionContext): DBIO[Unit] = {
//    if (useSlickUpserts) {
//      for {
//        _ <- dataAccess.summaryStatusEntryIdsAutoInc.
//          insertOrUpdate(SummaryStatusEntry(summaryName, newSummaryPosition))
//      } yield ()
//    } else {
//      for {
        dataAccess.summaryPositionForSummaryNameAndExpectedPosition((summaryName, oldSummaryPosition)).update(newSummaryPosition).map (updated =>

          if (updated == 0) {
            println(s"Summarizer did not move forwards to ${newSummaryPosition}. Expected old summary position was not met: ${oldSummaryPosition}")
            ()
          }
        )
//        _ <- updateCount match {
//          case 0 =>
//            dataAccess.summaryStatusEntryIdsAutoInc +=
//              SummaryStatusEntry(summaryName, newSummaryPosition)
//          case _ => assertUpdateCount("upsertSummaryStatusEntrySummaryPosition", updateCount, 1)
//        }
//      } yield ()
//    }
  }

  def getSummaryStatusIfNotBelowThreshold(summaryName: String,
                                          summaryPosition: Long)
                                         (implicit ec: ExecutionContext): DBIO[Option[Long]] = {
    for {
      value <- dataAccess.summaryPositionForSummaryNameIfAboveThreshold((summaryName, summaryPosition)).result
    } yield value.headOption
  }

  def ensureSummaryStatusIsBelowThreshold(summaryName: String,
                                          summaryPosition: Long)
                                         (implicit ec: ExecutionContext): DBIO[Boolean] = {
    for {
      update <- dataAccess.summaryPositionForSummaryNameIfAboveThreshold((summaryName, summaryPosition)).update(summaryPosition - 1)
    } yield update > 0
  }
}
