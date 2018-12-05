package centaur

import cats.effect.IO
import centaur.AbstractCromwellEngineOrBackendUpgradeTestCaseSpec._
import centaur.test.standard.CentaurTestCase
import cromwell.database.slick.{EngineSlickDatabase, MetadataSlickDatabase, SlickDatabase}
import cromwell.database.sql.SqlDatabase
import org.scalatest.{Assertions, BeforeAndAfter, DoNotDiscover}
import shapeless.syntax.typeable._

import scala.concurrent.Future

@DoNotDiscover
abstract class AbstractCromwellEngineOrBackendUpgradeTestCaseSpec(cromwellBackends: List[String])
  extends AbstractCentaurTestCaseSpec(cromwellBackends)
    with CentaurTestSuiteShutdown
    with BeforeAndAfter {

  def this() = this(CentaurTestSuite.cromwellBackends)

  private val cromwellDatabase = CromwellDatabase.fromConfig(CentaurConfig.conf)
  private val engineSlickDatabaseOption = cromwellDatabase.engineDatabase.cast[EngineSlickDatabase]
  private val metadataSlickDatabaseOption = cromwellDatabase.metadataDatabase.cast[MetadataSlickDatabase]

  override protected def beforeAll(): Unit = {
    super.beforeAll()
    val beforeAllIo = for {
      _ <- checkIsEmpty(cromwellDatabase.engineDatabase, cromwellDatabase.engineDatabase.existsJobKeyValueEntries(), testType)
      _ <- checkIsEmpty(cromwellDatabase.metadataDatabase, cromwellDatabase.metadataDatabase.existsMetadataEntries(), testType)
    } yield ()
    beforeAllIo.unsafeRunSync()
  }

  private def failNotSlick(database: SqlDatabase): IO[Unit] = {
    IO.raiseError(new RuntimeException(s"Expected a slick database for ${database.connectionDescription}."))
  }

  after {
    val afterIo = for {
      _ <- IO(CromwellManager.stopCromwell("Resetting Cromwell database back to prior version."))
      _ <- engineSlickDatabaseOption.fold(failNotSlick(cromwellDatabase.engineDatabase))(recreateDatabase)
      _ <- metadataSlickDatabaseOption.fold(failNotSlick(cromwellDatabase.metadataDatabase))(recreateDatabase)
      _ <- IO(CentaurTestSuite.startCromwell())
    } yield ()
    afterIo.unsafeRunSync()
  }

  allTestCases filter isMatchingUpgradeTest foreach executeStandardTest

  def testType: String

  def isMatchingUpgradeTest(testCase: CentaurTestCase): Boolean
}


object AbstractCromwellEngineOrBackendUpgradeTestCaseSpec {
  private def checkIsEmpty(database: SqlDatabase, lookup: => Future[Boolean], testType: => String): IO[Unit] = {
    IO.fromFuture(IO(lookup)).flatMap(exists =>
      if (exists) {
        IO(Assertions.fail(
          s"Database ${database.connectionDescription} contains data. " +
            s"$testType tests should only be run on a completely empty database. " +
            "You may need to manually drop and recreate the database to continue."
        ))
      } else {
        IO.unit
      }
    )
  }

  private def recreateDatabase(slickDatabase: SlickDatabase): IO[Unit] = {
    import slickDatabase.dataAccess.driver.api._
    val schemaName = slickDatabase.databaseConfig.getString("db.schema")
    //noinspection SqlDialectInspection
    for {
      _ <- IO.fromFuture(IO(slickDatabase.database.run(sqlu"""DROP SCHEMA IF EXISTS #$schemaName""")))
      _ <- IO.fromFuture(IO(slickDatabase.database.run(sqlu"""CREATE SCHEMA #$schemaName""")))
    } yield ()
  }
}
