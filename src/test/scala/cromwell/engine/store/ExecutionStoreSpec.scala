package cromwell.engine.store

import org.scalatest.{Matchers, FlatSpec}
import cromwell.binding._
import cromwell.engine.db._
import cromwell.engine.store.ExecutionStore.ExecutionStatus
import ExecutionStoreSpec._

object ExecutionStoreSpec {
  case class TestCallInfo(callFqn: FullyQualifiedName, status: CallStatus) extends CallInfo

  val Wdl =
    """
      |task a {
      |  command { echo '${message}' }
      |  output {
      |    String message = read_string(stdout())
      |    Int constant = 100
      |  }
      |}
      |task b {
      |  command { echo '${message} - ${Int integer}' }
      |}
      |workflow wf {
      |  call a
      |  call b {
      |    input: message=a.message, integer=a.constant - 75
      |  }
      |}
    """.stripMargin

  val Namespace = WdlNamespace.load(Wdl)

  val InitialStore: Set[CallInfo] = Set(TestCallInfo("a", ExecutionStatus.Done))
}

class ExecutionStoreSpec extends FlatSpec with Matchers {
  "An execution store" should "be startable in an arbitrary state" in {
    val store = ExecutionStore(Namespace, InitialStore)
    assert(!store.isWorkflowDone)
    val runnableCalls = store.startRunnableCalls
    runnableCalls should have size 1
    runnableCalls.head.taskFqn should be ("b")
  }
}
