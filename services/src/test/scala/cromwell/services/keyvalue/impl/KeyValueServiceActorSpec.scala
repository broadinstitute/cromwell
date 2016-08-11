package cromwell.services.keyvalue.impl

import akka.pattern._
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import cromwell.services.ServicesSpec
import cromwell.services.keyvalue.KeyValueServiceActor._

class KeyValueServiceActorSpec extends ServicesSpec("KeyValue") {

  val cromwellConfig = ConfigFactory.parseString(
    s"""services: {
        |  KeyValue: {
        |    class: "cromwell.services.keyvalue.KeyValueServiceActor"
        |    config {
        |      option1: "value1"
        |    }
        |  }
        |}
     """.stripMargin
  )

  val emptyConfig = ConfigFactory.empty()
  val sqlKvServiceActor = system.actorOf(SqlKeyValueServiceActor.props(emptyConfig, emptyConfig))
  val wfID = WorkflowId.randomId()

  val jobKey1 = KvJobKey("some_FQN", Option(-1), 1)
  val jobKey2 = KvJobKey("some_FQN", Option(-1), 2)

  val kvPair1 = KvPair(ScopedKey(wfID, jobKey1, "k1"), Option("v1"))
  val kvPair2 = KvPair(ScopedKey(wfID, jobKey1, "k2"), Option("v2"))
  val kvPair3 = KvPair(ScopedKey(wfID, jobKey2, "k1"), Option("v1"))

  "KeyValueServiceActor" should {
    "insert a key/value" in {
      val kvPut1 = KvPut(KvPair(ScopedKey(wfID, jobKey1, "k1"), Option("v1")))

          (for {
            putResult <- (sqlKvServiceActor ? kvPut1).mapTo[KvResponse]
            _ = putResult shouldEqual KvPutSuccess(kvPut1)
            putResult <- (sqlKvServiceActor ? KvPut(kvPair2)).mapTo[KvResponse]
            _ = putResult shouldEqual KvPutSuccess(KvPut(kvPair2))
            putResult <- (sqlKvServiceActor ? KvPut(kvPair3)).mapTo[KvResponse]
            _ = putResult shouldEqual KvPutSuccess(KvPut(kvPair3))
          } yield ()).futureValue
    }

    "return error if key doesn't exist" in {
      val scopedKeyNeverPut = ScopedKey(wfID, jobKey2, "k2")

          (for {
            getResult <- (sqlKvServiceActor ? KvGet(scopedKeyNeverPut)).mapTo[KvResponse]
            _ = getResult shouldEqual KvKeyLookupFailed(KvGet(scopedKeyNeverPut))
          } yield ()).futureValue
    }

    "be able to overwrite values" in {
      val kvPairOverride = KvPair(ScopedKey(wfID, jobKey1, "k1"), Option("v3"))

      (for {
        putResult <- (sqlKvServiceActor ? KvPut(kvPair2)).mapTo[KvResponse]
        _ = putResult shouldEqual KvPutSuccess(KvPut(kvPair2))
        putResult <- (sqlKvServiceActor ? KvPut(kvPairOverride)).mapTo[KvResponse]
        _ = putResult shouldEqual KvPutSuccess(KvPut(kvPairOverride))
      } yield ()).futureValue
    }
  }
}
