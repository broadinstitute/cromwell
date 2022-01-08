package cromwell.services.database

import cromwell.database.slick.SlickDatabase
import slick.jdbc.TransactionIsolation

import scala.concurrent.Future
import scala.concurrent.duration.Duration

/**
  * Exposes the protected method runTransaction via a runTestTransaction.
  */
trait TestSlickDatabase {
  slickDatabase: SlickDatabase =>

  import dataAccess.driver.api._

  def runTestTransaction[R](action: DBIO[R],
                            isolationLevel: TransactionIsolation = TransactionIsolation.RepeatableRead,
                            timeout: Duration = Duration.Inf,
                           ): Future[R] = {
    slickDatabase.runTransaction(action, isolationLevel, timeout)
  }
}
