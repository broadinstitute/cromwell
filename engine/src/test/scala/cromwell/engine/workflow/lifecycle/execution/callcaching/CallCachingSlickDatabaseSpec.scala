package cromwell.engine.workflow.lifecycle.execution.callcaching

import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.DbmsTest
import cromwell.core.WorkflowId
import cromwell.database.slick.EngineSlickDatabase
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables._
import cromwell.services.EngineServicesStore
import cromwell.services.ServicesStore.EnhancedSqlDatabase
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.ExecutionContext

class CallCachingSlickDatabaseSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {

  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(5, Seconds)), scaled(Span(100, Millis)))

  "SlickDatabase (hsqldb)" should behave like testWith("database")

  "SlickDatabase (mysql)" should behave like testWith("database-test-mysql")

  "SlickDatabase (mariadb)" should behave like testWith("database-test-mariadb")

  "SlickDatabase (postgresql)" should behave like testWith("database-test-postgresql")

  def testWith(configPath: String): Unit = {
    lazy val databaseConfig = ConfigFactory.load.getConfig(configPath)
    lazy val dataAccess = new EngineSlickDatabase(databaseConfig)
      .initialized(EngineServicesStore.EngineLiquibaseSettings)

    val idA = WorkflowId.randomId().toString
    val callA = "AwesomeWorkflow.GoodJob"
    val callCachingEntryA = CallCachingEntry(
      idA,
      callA,
      1,
      None,
      None,
      allowResultReuse = false
    )

    val callCachingHashEntriesA = Seq(
      CallCachingHashEntry(
        hashKey = "input: String s1",
        hashValue = "HASH_S1"
      ),
      CallCachingHashEntry(
        hashKey = "input: String s2",
        hashValue = "HASH_S2"
      ),
      CallCachingHashEntry(
        hashKey = "input: String s4",
        hashValue = "HASH_S4"
      )
    )

    val callCachingSimpletonsA = Seq(
      CallCachingSimpletonEntry("simpleKey", "simpleValue".toClobOption, "string")
    )

    val callCachingDetritusesA = Seq(
      CallCachingDetritusEntry("detritusKey", "detritusValue".toClobOption)
    )

    val aggregation = Option(CallCachingAggregationEntry("BASE_AGGREGATION", Option("FILE_AGGREGATION")))

    it should "honor allowResultReuse" taggedAs DbmsTest in {
      (for {
        _ <- dataAccess.addCallCaching(Seq(
          CallCachingJoin(
            callCachingEntryA,
            callCachingHashEntriesA,
            aggregation,
            callCachingSimpletonsA, callCachingDetritusesA
          )
        ),
          100
        )
        hasBaseAggregation <- dataAccess.hasMatchingCallCachingEntriesForBaseAggregation("BASE_AGGREGATION")
        _ = hasBaseAggregation shouldBe false
        hit <- dataAccess.findCacheHitForAggregation("BASE_AGGREGATION", Option("FILE_AGGREGATION"), callCachePathPrefixes = None, 1)
        _ = hit shouldBe empty
      } yield ()).futureValue
    }

    it should "retrieve CallCacheJoin for call" taggedAs DbmsTest in {
      (for {
        join <- dataAccess.callCacheJoinForCall(idA, callA, 1)
        _ = join shouldBe defined
        getJoin = join.get
        // We can't compare directly because the ones out from the DB have IDs filled in, so just compare the relevant values
        _ = getJoin
          .callCachingHashEntries
          .map(e => (e.hashKey, e.hashValue)) should contain theSameElementsAs callCachingHashEntriesA.map(e => (e.hashKey, e.hashValue))
        _ = getJoin
          .callCachingSimpletonEntries
          .map(e => (e.simpletonKey, e.simpletonValue.map(_.toRawString)))should contain theSameElementsAs callCachingSimpletonsA.map(e => (e.simpletonKey, e.simpletonValue.map(_.toRawString)))
        _ = getJoin
          .callCachingAggregationEntry
          .map(e => (e.baseAggregation, e.inputFilesAggregation)) shouldBe aggregation.map(e => (e.baseAggregation, e.inputFilesAggregation))
        _ = getJoin
          .callCachingDetritusEntries
          .map(e => (e.detritusKey, e.detritusValue.map(_.toRawString))) should contain theSameElementsAs callCachingDetritusesA.map(e => (e.detritusKey, e.detritusValue.map(_.toRawString)))
      } yield ()).futureValue
    }
  }
}
