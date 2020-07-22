package cromwell.engine.workflow.lifecycle.execution.stores

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus.{ExecutionStatus, _}
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.workflow.lifecycle.execution.stores.ExecutionStoreSpec._
import org.scalatest.{FlatSpec, Matchers}
import wom.graph._

class ExecutionStoreSpec extends FlatSpec with Matchers {



  it should "keep requiring an update until all 10000 call keys are started" in {

    def jobKeys: Map[JobKey, ExecutionStatus] = (0.until(10000).toList map {
      i => BackendJobDescriptorKey(noConnectionsGraphNode, Option(i), 1) -> NotStarted })
      .toMap

    var store: ExecutionStore = ActiveExecutionStore(jobKeys, needsUpdate = true)

    while (store.needsUpdate) {
      print("cycle")
      val update = store.update
      update.runnableKeys.size should be(1000)
      store = update.updatedStore.updateKeys(update.runnableKeys.map(_ -> QueuedInCromwell).toMap)
    }

    store.store(QueuedInCromwell).size should be(10000)
  }
}

object ExecutionStoreSpec {
  final case class MockJobKey(i: Int) extends CallKey {
    override def node: CallNode = noConnectionsGraphNode
    override def index: Option[Int] = Option(0)
    override def attempt: Int = 0
    override def tag: String = s"test_job_0"
  }

  val noConnectionsGraphNode: CommandCallNode = CommandCallNode(
    identifier = WomIdentifier("mock_task", "mock_wf.mock_task"),
    callable = null,
    inputPorts =  Set.empty[GraphNodePort.InputPort],
    inputDefinitionMappings = List.empty,
    nonInputBasedPrerequisites = Set.empty[GraphNode],
    outputIdentifierCompoundingFunction = (wi, _) => wi,
    sourceLocation = None
  )
}
