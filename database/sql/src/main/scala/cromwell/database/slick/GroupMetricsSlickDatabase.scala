package cromwell.database.slick

import cromwell.database.sql.GroupMetricsSqlDatabase
import cromwell.database.sql.tables.GroupMetricsEntry
import org.postgresql.util.PSQLException
import slick.jdbc.TransactionIsolation

import java.sql.{SQLIntegrityConstraintViolationException, Timestamp}
import scala.concurrent.{ExecutionContext, Future}

trait GroupMetricsSlickDatabase extends GroupMetricsSqlDatabase {
  this: EngineSlickDatabase =>

  import dataAccess.driver.api._

//  private val mutex = new Object

  override def recordGroupMetricsEntry(
    groupMetricsEntry: GroupMetricsEntry
  )(implicit ec: ExecutionContext): Future[Unit] = {

    // original code
//    val action = for {
//      updateCount <- dataAccess
//        .quotaExhaustionForGroupId(groupMetricsEntry.groupId)
//        .update(groupMetricsEntry.quotaExhaustionDetected)
//      _ <- updateCount match {
//        case 0 => dataAccess.groupMetricsEntryIdsAutoInc += groupMetricsEntry
//        case _ => assertUpdateCount("recordGroupMetricsEntry", updateCount, 1)
//      }
//    } yield ()

//    val action = for {
//      groupMetricsEntries <- dataAccess
//        .entryForGroupIdForUpdate(groupMetricsEntry.groupId)
//        .result
//      _ <- groupMetricsEntries.map {entry =>
//        for {
//          updateCount <- dataAccess.quotaExhaustionForGroupId(entry.groupId).update(groupMetricsEntry.quotaExhaustionDetected)
//          _ <- assertUpdateCount("recordGroupMetricsEntry", updateCount, 1)
//        } yield ()
//      }
//
//    } yield ()

//    mutex.synchronized {

    // working forUpdate but still deadlocks
//    val action = for {
//      _ <- dataAccess.entryForGroupIdForUpdate(groupMetricsEntry.groupId).result
//      updateCount <- dataAccess.quotaExhaustionForGroupId(groupMetricsEntry.groupId)
//        .update(groupMetricsEntry.quotaExhaustionDetected)
//      _ <- updateCount match {
//        case 0 => dataAccess.groupMetricsEntryIdsAutoInc += groupMetricsEntry
//        case _ => assertUpdateCount("recordGroupMetricsEntry", updateCount, 1)
//      }
//    } yield ()

    // uses upsert
//    val action = if(useSlickUpserts) {
//      for {
//        existingRow <- dataAccess.entryForGroupIdForUpdate(groupMetricsEntry.groupId).result.headOption
//        row = existingRow.map(_.copy(quotaExhaustionDetected = groupMetricsEntry.quotaExhaustionDetected)) getOrElse GroupMetricsEntry(groupMetricsEntry.groupId, groupMetricsEntry.quotaExhaustionDetected)
//        _ <- dataAccess.groupMetricsEntryIdsAutoInc.insertOrUpdate(row)
//      } yield ()
//    } else {
//      for {
//              updateCount <- dataAccess
//                .quotaExhaustionForGroupId(groupMetricsEntry.groupId)
//                .update(groupMetricsEntry.quotaExhaustionDetected)
//              _ <- updateCount match {
//                case 0 => dataAccess.groupMetricsEntryIdsAutoInc += groupMetricsEntry
//                case _ => assertUpdateCount("recordGroupMetricsEntry", updateCount, 1)
//              }
//            } yield ()
//    }

    // insert first; if it fails try update
    val action = for {
      _ <- dataAccess.groupMetricsEntryIdsAutoInc += groupMetricsEntry
    } yield ()

    runTransaction(action, TransactionIsolation.ReadCommitted).flatMap { _ =>
      Future.successful(())
    }.recoverWith {
      case _: SQLIntegrityConstraintViolationException =>
        println(s"#### FIND ME: matched with SQLIntegrityConstraintViolationException")
        val updateAction = for {
          _ <- dataAccess
            .quotaExhaustionForGroupId(groupMetricsEntry.groupId)
            .update(groupMetricsEntry.quotaExhaustionDetected)
        } yield ()
        runTransaction(updateAction, TransactionIsolation.ReadCommitted)

      case e: PSQLException if e.getMessage.contains("duplicate key value violates unique constraint \"UC_GROUP_METRICS_ENTRY_GI\"") =>
        println(s"#### FIND ME: matched with PSQLException")
        val updateAction = for {
          _ <- dataAccess
            .quotaExhaustionForGroupId(groupMetricsEntry.groupId)
            .update(groupMetricsEntry.quotaExhaustionDetected)
        } yield ()
        runTransaction(updateAction, TransactionIsolation.ReadCommitted)

      case e =>  Future.failed(e)
    }




//      runTransaction(action, TransactionIsolation.ReadCommitted) // handles deadlock but doesn't keep group ID unique; inserts multiple rows for same group id in parallel case
//    runTransaction(action, TransactionIsolation.Serializable) // doesn't handle deadlock in all MySql and Postgres versions
//    runTransaction(action, TransactionIsolation.RepeatableRead) // doesn't handle deadlock in MySQL 5.6, and both Postgres versions (MySql 8 works fine)
//    }
  }

  override def countGroupMetricsEntries(groupId: String)(implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.countGroupMetricsEntriesForGroupId(groupId).result
    runTransaction(action)
  }

  override def getQuotaExhaustedGroups(
    thresholdTimestamp: Timestamp
  )(implicit ec: ExecutionContext): Future[Seq[String]] = {
    val action = dataAccess.groupsExperiencingQuotaExhaustion(thresholdTimestamp).result
    runTransaction(action)
  }
}
