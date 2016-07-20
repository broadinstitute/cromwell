package cromwell.engine

import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore._
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.Matchers

import scala.concurrent.duration._
import scala.language.postfixOps
import scalaz.NonEmptyList

class WorkflowStoreActorSpec extends CromwellTestkitSpec with Matchers {
  val helloWorldSourceFiles = HelloWorld.asWorkflowSources()

  /**
    * Fold down a list of WorkflowToStart's, checking that their IDs are all unique
    */
  private def checkDistinctIds(list: Traversable[WorkflowToStart]): Boolean = {
    def folderFunction(knownDistinct: (List[WorkflowToStart], Boolean), next: WorkflowToStart) = {
      val (list, distinct) = knownDistinct
      if (!distinct) {
        (list :+ next, false)
      }
      else {
        (list :+ next, !list.map(_.id).contains(next.id))
      }
    }

    list.foldLeft((List.empty[WorkflowToStart], true))(folderFunction)._2
  }

  "The WorkflowStoreActor" should {
    "return an ID for a submitted workflow" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store))
      storeActor ! SubmitWorkflow(helloWorldSourceFiles)
      expectMsgType[WorkflowSubmittedToStore](10 seconds)
    }

    "return 3 IDs for a batch submission of 3" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store))
      storeActor ! BatchSubmitWorkflows(NonEmptyList(helloWorldSourceFiles, helloWorldSourceFiles, helloWorldSourceFiles))
      expectMsgPF(10 seconds) {
        case WorkflowsBatchSubmittedToStore(ids) => ids.size shouldBe 3
      }
    }

    "fetch exactly N workflows" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store))
      storeActor ! BatchSubmitWorkflows(NonEmptyList(helloWorldSourceFiles, helloWorldSourceFiles, helloWorldSourceFiles))
      val insertedIds = expectMsgType[WorkflowsBatchSubmittedToStore](10 seconds).workflowIds.list


      storeActor ! FetchRunnableWorkflows(2)
      expectMsgPF(10 seconds) {
        case NewWorkflowsToStart(workflowNel) =>
          workflowNel.size shouldBe 2
          checkDistinctIds(workflowNel.list) shouldBe true
          workflowNel.foreach {
            case WorkflowToStart(id, sources, state) =>
              insertedIds.contains(id) shouldBe true
              sources shouldBe helloWorldSourceFiles
              state shouldBe Submitted
          }
      }
    }

    "return only the remaining workflows if N is larger than size" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store))
      storeActor ! BatchSubmitWorkflows(NonEmptyList(helloWorldSourceFiles, helloWorldSourceFiles, helloWorldSourceFiles))
      val insertedIds = expectMsgType[WorkflowsBatchSubmittedToStore](10 seconds).workflowIds.list


      storeActor ! FetchRunnableWorkflows(100)
      expectMsgPF(10 seconds) {
        case NewWorkflowsToStart(workflowNel) =>
          workflowNel.size shouldBe 3
          checkDistinctIds(workflowNel.list) shouldBe true
          workflowNel.foreach {
            case WorkflowToStart(id, sources, state) =>
              insertedIds.contains(id) shouldBe true
              sources shouldBe helloWorldSourceFiles
              state shouldBe Submitted
          }
      }
    }

    "remove workflows which exist" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store))
      storeActor ! SubmitWorkflow(helloWorldSourceFiles)
      val id = expectMsgType[WorkflowSubmittedToStore](10 seconds).workflowId
      storeActor ! RemoveWorkflow(id)
      storeActor ! FetchRunnableWorkflows(100)
      expectMsgPF(10 seconds) {
        case NoNewWorkflowsToStart => // Great
        case x => fail(s"Unexpected response from supposedly empty WorkflowStore: $x")
      }
    }

    "remain responsive if you ask to remove a workflow it doesn't have" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store))
      val id = WorkflowId.randomId()
      storeActor ! RemoveWorkflow(id)

      storeActor ! FetchRunnableWorkflows(100)
      expectMsgPF(10 seconds) {
        case NoNewWorkflowsToStart => // Great
        case x => fail(s"Unexpected response from supposedly empty WorkflowStore: $x")
      }
    }
  }
}
