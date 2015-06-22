package cromwell.engine.store

import cromwell.util.SampleWdl.SubtractionWorkflow
import org.scalatest.{Matchers, FlatSpec}
import cromwell.binding._
import cromwell.engine.db._
import cromwell.engine.store.ExecutionStore.ExecutionStatus
import ExecutionStoreSpec._

object ExecutionStoreSpec {
  case class TestCallInfo(callFqn: FullyQualifiedName, status: CallStatus) extends CallInfo
  val Namespace = WdlNamespace.load(SubtractionWorkflow.WdlSource)
  val InitialStore: Set[CallInfo] = Set(TestCallInfo("a", ExecutionStatus.Done))
}

class ExecutionStoreSpec extends FlatSpec with Matchers {
  "An execution store" should "be startable in an arbitrary state" in {
    val store = ExecutionStore(Namespace, InitialStore)
    store.isWorkflowDone should be (false)
    val runnableCalls = store.startRunnableCalls
    runnableCalls should have size 1
    runnableCalls.head.taskFqn should be ("b")
  }
}
