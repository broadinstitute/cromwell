package cromwell.engine.workflow.lifecycle.execution.stores

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus.{ExecutionStatus, _}
import cromwell.core.JobKey
import cromwell.engine.workflow.lifecycle.execution.stores.ExecutionStoreSpec._
import org.scalatest.{BeforeAndAfter, FlatSpec, Matchers}
import wom.graph._

import scala.util.Random

class ExecutionStoreSpec extends FlatSpec with Matchers with BeforeAndAfter {

  var store: ExecutionStore = _

  before {
    def jobKeys: Map[JobKey, ExecutionStatus] = (0.until(10000).toList map {
      i => BackendJobDescriptorKey(noConnectionsGraphNode, Option(i), 1) -> NotStarted })
      .toMap

    store = ActiveExecutionStore(jobKeys, needsUpdate = true)
  }

  def updateStoreToEnqueueNewlyRunnableJobs(): Unit = {
    while (store.needsUpdate) {
      val update = store.update
      store = update.updatedStore.updateKeys(update.runnableKeys.map(_ -> QueuedInCromwell).toMap)
    }
  }

  it should "keep allowing updates while 10000 call keys are enqueued and then started in small batches" in {

    updateStoreToEnqueueNewlyRunnableJobs()

    store.store.getOrElse(QueuedInCromwell, List.empty).size should be(10000)
    store.store.getOrElse(Running, List.empty).size should be(0)

    while(store.store.getOrElse(Running, List.empty).size < 10000) {
      val newlyRunning = Random.nextInt(1000)
      store = store.updateKeys(store.store(QueuedInCromwell).take(newlyRunning).map(j => j -> Running).toMap)
      updateStoreToEnqueueNewlyRunnableJobs()
    }

    store.store.getOrElse(QueuedInCromwell, List.empty).size should be(0)
    store.store.getOrElse(Running, List.empty).size should be(10000)
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
