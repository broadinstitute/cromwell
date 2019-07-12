package cromwell.services.database

import better.files._
import cromwell.core.Tags._
import org.scalactic.StringNormalizations
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

/**
  * Ensures the in-memory HSQLDB is setup with MVCC transaction isolation.
  *
  * http://www.hsqldb.org/doc/guide/sessions-chapt.html#snc_tx_mvcc
  * https://en.wikipedia.org/wiki/Multiversion_concurrency_control
  */
class HsqldbTransactionIsolationSpec extends FlatSpec with Matchers with ScalaFutures with StringNormalizations {

  CromwellDatabaseType.All foreach { databaseType =>
    behavior of s"HSQLDB transaction isolation for ${databaseType.name}"

    it should "be mvcc" taggedAs DbmsTest in {
      for {
        slickDatabase <- DatabaseTestKit.initializedDatabaseFromSystem(databaseType, HsqldbDatabaseSystem).autoClosed
      } {
        import slickDatabase.dataAccess.driver.api._
        //noinspection SqlDialectInspection
        val getHsqldbTx =
          sql"""SELECT PROPERTY_VALUE
                FROM INFORMATION_SCHEMA.SYSTEM_PROPERTIES
                WHERE PROPERTY_NAME = 'hsqldb.tx'
             """.as[String].head

        val future = slickDatabase.database.run(getHsqldbTx)
        (future.futureValue shouldEqual "mvcc") (after being lowerCased)
      }
    }
  }
}
