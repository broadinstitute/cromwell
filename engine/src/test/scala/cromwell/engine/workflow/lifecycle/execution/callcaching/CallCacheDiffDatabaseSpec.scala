package cromwell.engine.workflow.lifecycle.execution.callcaching

import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.DbmsTest
import cromwell.core.WorkflowId
import cromwell.database.slick.SlickDatabase
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.{CallCachingEntry, CallCachingHashEntry}
import cromwell.services.ServicesStore
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.specs2.mock.Mockito

import scala.concurrent.ExecutionContext

class CallCacheDiffDatabaseSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {
  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(5, Seconds)), scaled(Span(100, Millis)))
  
  "SlickDatabase (hsqldb)" should behave like testWith("database")

  "SlickDatabase (mysql)" should behave like testWith("database-test-mysql")

  def testWith(configPath: String): Unit = {
    import ServicesStore.EnhancedSqlDatabase

    lazy val databaseConfig = ConfigFactory.load.getConfig(configPath)
    lazy val dataAccess = new SlickDatabase(databaseConfig).initialized

    val callCachingEntryA = CallCachingEntry(
      WorkflowId.randomId().toString,
      "AwesomeWorkflow.GoodJob",
      1,
      None,
      None,
      allowResultReuse = true
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
    
    val callCachingEntryB = CallCachingEntry(
      WorkflowId.randomId().toString,
      "BetterWorkflow.GreatJob",
      1,
      None,
      None,
      allowResultReuse = true
    )

    val callCachingHashEntriesB = Seq(
      CallCachingHashEntry(
        hashKey = "input: String s1",
        hashValue = "HASH_S1"
      ),
      CallCachingHashEntry(
        hashKey = "input: String s2",
        hashValue = "HASH_NOT_S2"
      ),
      CallCachingHashEntry(
        hashKey = "input: String s3",
        hashValue = "HASH_S3"
      )
    )
    
    val callCachingEntryC = CallCachingEntry(
      WorkflowId.randomId().toString,
      "SameWorkflow.SameJob",
      1,
      None,
      None,
      allowResultReuse = true
    )

    val callCachingHashEntriesC = callCachingHashEntriesB

    it should "find correct hash diff" taggedAs DbmsTest in {
      (for {
        _ <- dataAccess.addCallCaching(Seq(
          CallCachingJoin(
            callCachingEntryA,
            callCachingHashEntriesA,
            None, Seq.empty, Seq.empty
          ),
          CallCachingJoin(
            callCachingEntryB,
            callCachingHashEntriesB,
            None, Seq.empty, Seq.empty
          )
        ),
          100
        )
        hashDiff <- dataAccess.diffCallCacheHashes(
          callCachingEntryA.workflowExecutionUuid, callCachingEntryA.callFullyQualifiedName, callCachingEntryA.jobIndex,
          callCachingEntryB.workflowExecutionUuid, callCachingEntryB.callFullyQualifiedName, callCachingEntryB.jobIndex
        )
        _ = hashDiff.cacheEntryA.workflowExecutionUuid shouldBe callCachingEntryA.workflowExecutionUuid
        _ = hashDiff.cacheEntryA.callFullyQualifiedName shouldBe callCachingEntryA.callFullyQualifiedName
        _ = hashDiff.cacheEntryA.jobIndex shouldBe callCachingEntryA.jobIndex
        _ = hashDiff.cacheEntryA.allowResultReuse shouldBe callCachingEntryA.allowResultReuse
      
        _ = hashDiff.cacheEntryB.workflowExecutionUuid shouldBe callCachingEntryB.workflowExecutionUuid
        _ = hashDiff.cacheEntryB.callFullyQualifiedName shouldBe callCachingEntryB.callFullyQualifiedName
        _ = hashDiff.cacheEntryB.jobIndex shouldBe callCachingEntryB.jobIndex
        _ = hashDiff.cacheEntryB.allowResultReuse shouldBe callCachingEntryB.allowResultReuse
      
        _ = hashDiff.diff should contain theSameElementsAs List(
          Option("input: String s2" -> "HASH_S2") -> Option("input: String s2" -> "HASH_NOT_S2"),
          Option("input: String s4" -> "HASH_S4") -> None,
          None -> Option("input: String s3" -> "HASH_S3")
        )
      } yield ()).futureValue
    }

    it should "return empty diff if calls have same hashes" taggedAs DbmsTest in {
      (for {
        _ <- dataAccess.addCallCaching(Seq(
          CallCachingJoin(
            callCachingEntryC,
            callCachingHashEntriesC,
            None, Seq.empty, Seq.empty
          )
        ),
          100
        )
        hashDiff <- dataAccess.diffCallCacheHashes(
          callCachingEntryB.workflowExecutionUuid, callCachingEntryB.callFullyQualifiedName, callCachingEntryB.jobIndex,
          callCachingEntryC.workflowExecutionUuid, callCachingEntryC.callFullyQualifiedName, callCachingEntryC.jobIndex
        )
        _ = hashDiff.cacheEntryA.workflowExecutionUuid shouldBe callCachingEntryB.workflowExecutionUuid
        _ = hashDiff.cacheEntryA.callFullyQualifiedName shouldBe callCachingEntryB.callFullyQualifiedName
        _ = hashDiff.cacheEntryA.jobIndex shouldBe callCachingEntryB.jobIndex
        _ = hashDiff.cacheEntryA.allowResultReuse shouldBe callCachingEntryB.allowResultReuse

        _ = hashDiff.cacheEntryB.workflowExecutionUuid shouldBe callCachingEntryC.workflowExecutionUuid
        _ = hashDiff.cacheEntryB.callFullyQualifiedName shouldBe callCachingEntryC.callFullyQualifiedName
        _ = hashDiff.cacheEntryB.jobIndex shouldBe callCachingEntryC.jobIndex
        _ = hashDiff.cacheEntryB.allowResultReuse shouldBe callCachingEntryC.allowResultReuse

        _ = hashDiff.diff shouldBe empty
      } yield ()).futureValue
    }

    it should "return empty diff for same call" taggedAs DbmsTest in {
      (for {
        hashDiff <- dataAccess.diffCallCacheHashes(
          callCachingEntryA.workflowExecutionUuid, callCachingEntryA.callFullyQualifiedName, callCachingEntryA.jobIndex,
          callCachingEntryA.workflowExecutionUuid, callCachingEntryA.callFullyQualifiedName, callCachingEntryA.jobIndex
        )
        _ = hashDiff.cacheEntryA.workflowExecutionUuid shouldBe callCachingEntryA.workflowExecutionUuid
        _ = hashDiff.cacheEntryA.callFullyQualifiedName shouldBe callCachingEntryA.callFullyQualifiedName
        _ = hashDiff.cacheEntryA.jobIndex shouldBe callCachingEntryA.jobIndex
        _ = hashDiff.cacheEntryA.allowResultReuse shouldBe callCachingEntryA.allowResultReuse

        _ = hashDiff.cacheEntryB.workflowExecutionUuid shouldBe callCachingEntryA.workflowExecutionUuid
        _ = hashDiff.cacheEntryB.callFullyQualifiedName shouldBe callCachingEntryA.callFullyQualifiedName
        _ = hashDiff.cacheEntryB.jobIndex shouldBe callCachingEntryA.jobIndex
        _ = hashDiff.cacheEntryB.allowResultReuse shouldBe callCachingEntryA.allowResultReuse

        _ = hashDiff.diff shouldBe empty
      } yield ()).futureValue
    }

    it should "fail properly if the cache entry for call A is not found" taggedAs DbmsTest in {
      val shouldFail = for {
        _ <- dataAccess.diffCallCacheHashes(
          "does not", "exist", 0,
          callCachingEntryB.workflowExecutionUuid, callCachingEntryB.callFullyQualifiedName, callCachingEntryB.jobIndex
        )
      } yield ()
      
      ScalaFutures.whenReady(shouldFail.failed) { e =>
        e shouldBe an [Exception]
        e.getMessage shouldBe "Cannot find a cache entry for does not:exist:0"
      }
    }

    it should "fail properly if the cache entry for call B is not found" taggedAs DbmsTest in {
      val shouldFail = for {
        _ <- dataAccess.diffCallCacheHashes(
          callCachingEntryB.workflowExecutionUuid, callCachingEntryB.callFullyQualifiedName, callCachingEntryB.jobIndex,
          "does not", "existB", 0
        )
      } yield ()

      ScalaFutures.whenReady(shouldFail.failed) { e =>
        e shouldBe an [Exception]
        e.getMessage shouldBe "Cannot find a cache entry for does not:existB:0"
      }
    }

    it should "fail properly if none of the cache entries are found " taggedAs DbmsTest in {
      val shouldFail = for {
        _ <- dataAccess.diffCallCacheHashes(
          "does not", "exist", 0,
          "does not", "exist either", 0
        )
      } yield ()

      ScalaFutures.whenReady(shouldFail.failed) { e =>
        e shouldBe an [Exception]
        e.getMessage shouldBe "Cannot find cache entries for does not:exist:0, does not:exist either:0"
      }
    }
  }
}
