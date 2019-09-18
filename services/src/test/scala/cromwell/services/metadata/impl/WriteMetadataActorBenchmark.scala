package cromwell.services.metadata.impl

import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.core.Tags.IntegrationTest
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.database.slick.MetadataSlickDatabase
import cromwell.services.MetadataServicesStore
import cromwell.services.database.{DatabaseTestKit, EngineDatabaseType, MetadataDatabaseType, MysqlEarliestDatabaseSystem}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
class WriteMetadataActorBenchmark extends TestKitSuite with FlatSpecLike with Eventually with Matchers {
  override implicit val patienceConfig = PatienceConfig(scaled(30.seconds), 1.second)

  behavior of "WriteMetadataActor"
  
  val workflowId = WorkflowId.randomId()
  val registry = TestProbe().ref

  def makeEvent = {
    MetadataEvent(MetadataKey(workflowId, None, "metadata_key"), MetadataValue("My Value"))
  }

  def time[T](description: String)(thunk: => T): T = {
    val t1 = System.currentTimeMillis
    val x = thunk
    val t2 = System.currentTimeMillis
    println(description + " took " + (t2 - t1) + " ms")
    x
  }

  private lazy val dataAccess = new MetadataDatabaseAccess with MetadataServicesStore {
    override val metadataDatabaseInterface: MetadataSlickDatabase = {
      val databaseSystem = MysqlEarliestDatabaseSystem
      // NOTE: EngineLiquibaseSettings **MUST** always run before the MetadataLiquibaseSettings
      DatabaseTestKit.initializedDatabaseFromSystem(EngineDatabaseType, databaseSystem).close()
      DatabaseTestKit.initializedDatabaseFromSystem(MetadataDatabaseType, databaseSystem)
    }
  }

  it should "provide good throughput" taggedAs IntegrationTest in {
    val writeActor = TestFSMRef(new WriteMetadataActor(1000, 5.seconds, registry, Int.MaxValue) {
      override val metadataDatabaseInterface: MetadataSlickDatabase = dataAccess.metadataDatabaseInterface
    })
    
    time("metadata write") {
      (0 to 1 * 1000 * 1000)
        .map(_ => makeEvent)
        .grouped(100)
        .map(PutMetadataAction(_))
        .foreach(writeActor.!)

      eventually {
        writeActor.underlyingActor.stateData.weight shouldBe 0
      }
    }
  }

  it should "close the database" taggedAs IntegrationTest in {
    dataAccess.metadataDatabaseInterface.close()
  }
}
