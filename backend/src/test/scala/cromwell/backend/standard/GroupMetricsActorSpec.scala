package cromwell.backend.standard

import akka.actor.ActorSystem
import akka.testkit.{TestActorRef, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.standard.GroupMetricsActor.{
  GetQuotaExhaustedGroups,
  GetQuotaExhaustedGroupsSuccess,
  LogQuotaExhaustedGroups,
  RecordGroupQuotaExhaustion
}
import cromwell.database.slick.EngineSlickDatabase
import cromwell.database.sql.tables.GroupMetricsEntry
import cromwell.services.EngineServicesStore
import cromwell.services.ServicesStore.EnhancedSqlDatabase
import org.scalatest.concurrent.Eventually.eventually
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.sql.Timestamp
import scala.concurrent.duration.DurationInt
import scala.concurrent.{ExecutionContext, Future}

class GroupMetricsActorSpec extends AnyFlatSpec with Matchers {

  implicit val system: ActorSystem = ActorSystem("GroupMetricsActorSpec")

  val testHogGroup: String = "groot-hog-group"
  val DatabaseConfig: Config = ConfigFactory.load.getConfig("database")
  var recordMethodCallCount: Int = 0

  def databaseInterface(): EngineSlickDatabase =
    new EngineSlickDatabase(DatabaseConfig) {
      override def recordGroupMetricsEntry(groupMetricsEntry: GroupMetricsEntry)(implicit
        ec: ExecutionContext
      ): Future[Unit] = {
        recordMethodCallCount = recordMethodCallCount + 1
        Future.successful(())
      }

      override def getQuotaExhaustedGroups(thresholdTimestamp: Timestamp)(implicit
        ec: ExecutionContext
      ): Future[Seq[String]] =
        Future.successful(List(testHogGroup))
    }.initialized(EngineServicesStore.EngineLiquibaseSettings)

  behavior of "GroupMetricsActor"

  it should "receive new quota exhaustion message and call database function" in {
    val db = databaseInterface()
    val mockGroupMetricsActor = TestActorRef(GroupMetricsActor.props(db, 15, 5.minutes))

    mockGroupMetricsActor.tell(RecordGroupQuotaExhaustion(testHogGroup), TestProbe().ref)

    eventually {
      recordMethodCallCount shouldBe 1
    }
  }

  it should "respond with groups in quota exhaustion" in {
    val db = databaseInterface()
    val mockGroupMetricsActor = TestActorRef(GroupMetricsActor.props(db, 15, 5.minutes))
    val requestActor = TestProbe()

    mockGroupMetricsActor.tell(GetQuotaExhaustedGroups, requestActor.ref)

    requestActor.expectMsgPF(5.seconds) { case GetQuotaExhaustedGroupsSuccess(expectedList) =>
      expectedList shouldBe List(testHogGroup)
    }
  }
}
