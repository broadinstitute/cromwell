package cromwell.engine.workflow.lifecycle.execution.stores

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.ExecutionStatus.{ExecutionStatus, _}
import cromwell.core.JobKey
import cromwell.engine.workflow.lifecycle.execution.stores.ExecutionStoreSpec._
import org.scalatest.BeforeAndAfter
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wom.graph._

import scala.util.Random

class ExecutionStoreSpec extends AnyFlatSpec with Matchers with BeforeAndAfter {

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

  it should "keep allowing updates while 10000 call keys are enqueued and then started" in {

    updateStoreToEnqueueNewlyRunnableJobs()

    store.store(QueuedInCromwell).size should be(1000)
    store.store(WaitingForQueueSpace).size should be(9000)
    store.store.contains(Running) should be(false)

    var iteration = 0
    while(store.store.getOrElse(Running, List.empty).size < 10000) {
      store = store.updateKeys(store.store(QueuedInCromwell).map(j => j -> Running).toMap)
      updateStoreToEnqueueNewlyRunnableJobs()

      // In all situations, the queue should never go above the queue limit:
      (store.store.getOrElse(QueuedInCromwell, List.empty).size <= 1000) should be(true)

      // Additionally, as long as elements are waiting for queue space, the queue should be exactly 1000 long
      if (store.store.contains(WaitingForQueueSpace)) { store.store(QueuedInCromwell).size should be(1000) }

      iteration = iteration + 1
      store.store.getOrElse(Running, List.empty).size should be(iteration * 1000)
    }
  }

  it should "never enqueue more than the total allowed queue space" in {

    updateStoreToEnqueueNewlyRunnableJobs()

    store.store(QueuedInCromwell).size should be(1000)
    store.store(WaitingForQueueSpace).size should be(9000)
    store.store.contains(Running) should be(false)
    var currentlyRunning = 0

    while(store.store.getOrElse(Running, List.empty).size < 10000) {
      val newlyRunning: Iterable[JobKey] = store.store(QueuedInCromwell).take(Random.nextInt(1000))
      store = store.updateKeys(newlyRunning.map(j => j -> Running).toMap)
      updateStoreToEnqueueNewlyRunnableJobs()

      // In all situations, the queue should never go above the queue limit:
      (store.store.getOrElse(QueuedInCromwell, List.empty).size <= 1000) should be(true)

      // Additionally, as long as elements are waiting for queue space, the queue should be exactly 1000 long
      if (store.store.contains(WaitingForQueueSpace)) { store.store(QueuedInCromwell).size should be(1000) }

      store.store.getOrElse(Running, List.empty).size should be(currentlyRunning + newlyRunning.size)
      currentlyRunning = currentlyRunning + newlyRunning.size
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
