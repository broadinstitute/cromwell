package cromwell.backend.standard

import akka.actor.ActorSystem
import akka.testkit.{TestActorRef, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.standard.GroupMetricsActor.RecordGroupQuotaExhaustion
import cromwell.database.slick.EngineSlickDatabase
import cromwell.database.sql.tables.GroupMetricsEntry
import cromwell.services.EngineServicesStore
import cromwell.services.ServicesStore.EnhancedSqlDatabase
import org.scalatest.concurrent.Eventually.eventually
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

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
    }.initialized(EngineServicesStore.EngineLiquibaseSettings)

  behavior of "GroupMetricsActor"

  it should "receive new quota exhaustion message and call database function" in {
    val db = databaseInterface()
    val mockGroupMetricsActor = TestActorRef(GroupMetricsActor.props(db))

    mockGroupMetricsActor.tell(RecordGroupQuotaExhaustion(testHogGroup), TestProbe().ref)

    eventually {
      recordMethodCallCount shouldBe 1
    }
  }
}
