package cromwell.engine.workflow.lifecycle.execution.stores

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus.{ExecutionStatus, _}
import cromwell.core.JobKey
import cromwell.engine.workflow.lifecycle.execution.stores.ExecutionStoreSpec._
import org.scalatest.{FlatSpec, Matchers}
import wom.graph._

class ExecutionStoreSpec extends FlatSpec with Matchers {



  it should "keep requiring an update until all 10000 call keys are started" in {

    def jobKeys: Map[JobKey, ExecutionStatus] = (0.until(10000).toList map {
      i => BackendJobDescriptorKey(noConnectionsGraphNode, Option(i), 1) -> NotStarted })
      .toMap

    var store: ExecutionStore = ActiveExecutionStore(jobKeys, needsUpdate = true)

    def updateStoreToEnqueueNewlyRunnableJobs(): Unit = {
      while (store.needsUpdate) {
        val update = store.update
        store = update.updatedStore.updateKeys(update.runnableKeys.map(_ -> QueuedInCromwell).toMap)
      }
    }

    updateStoreToEnqueueNewlyRunnableJobs()

    store.store(QueuedInCromwell).size should be(1000)
    store.store(WaitingForQueueSpace).size should be(9000)

    while(store.store.getOrElse(Running, List.empty).size < 10000) {
      store = store.updateKeys(store.store(QueuedInCromwell).map(j => j -> Running).toMap)
      updateStoreToEnqueueNewlyRunnableJobs()
    }
  }
}

object ExecutionStoreSpec {

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
