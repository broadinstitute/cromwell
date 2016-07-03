package cromwell.engine

import cromwell.core.WorkflowId
import cromwell.engine.WorkflowStoreSpec.SpecWorkflowStore
import cromwell.engine.workflow.WorkflowStore
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.{FlatSpec, Matchers}

import scalaz.NonEmptyList

object WorkflowStoreSpec {
  class SpecWorkflowStore extends WorkflowStore
}

class WorkflowStoreSpec extends FlatSpec with Matchers {
  val sources = HelloWorld.asWorkflowSources()

  "A WorkflowStore" should "return an ID for a submitted workflow" in {
    val store = new SpecWorkflowStore
    val ids = store.add(NonEmptyList(sources))
    val storedIds = store.workflowStore map { _.id }
    storedIds should have size 1
    (storedIds diff ids.list) shouldBe empty
  }

  it should "return 3 IDs for a batch submission of 3" in {
    val store = new SpecWorkflowStore
    val newIds = store.add(NonEmptyList(sources, sources, sources))
    newIds should have size 3
    val storedIds = store.workflowStore map { _.id }
    storedIds should have size 3
    (storedIds diff newIds.list) shouldBe empty
  }

  it should "fetch no more than N workflows" in {
    val store = new SpecWorkflowStore
    val insertedIds = store.add(NonEmptyList(sources, sources, sources))
    val newIds = store.fetchRunnableWorkflows(1) map { _.id }
    newIds should have size 1
    val currentIds = store.workflowStore filter { _.state.isStartable } map { _.id }
    currentIds should be (insertedIds.list diff newIds)
  }

  it should "return only the remaining workflows if N is larger than size" in {
    val store = new SpecWorkflowStore
    val insertedIds = store.add(NonEmptyList(sources, sources, sources))
    val newIds = store.fetchRunnableWorkflows(store.workflowStore.size + 5) map { _.id }
    newIds should have size insertedIds.size
    val currentIds = store.workflowStore filter { _.state.isStartable } map { _.id }
    currentIds shouldBe empty
  }

  it should "remove workflows which exist" in {
    val store = new SpecWorkflowStore
    val id = store.add(NonEmptyList(sources)).head
    store.remove(id) shouldBe true
    store.workflowStore shouldBe empty
  }

  it should "not be super happy if you ask to remove a workflow it doesn't have" in {
    val store = new SpecWorkflowStore
    val id = store.add(NonEmptyList(sources)).head
    store.remove(WorkflowId.randomId()) shouldBe false
    store.workflowStore should have size 1
  }
}
