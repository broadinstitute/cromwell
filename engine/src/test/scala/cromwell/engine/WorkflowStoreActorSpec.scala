package cromwell.engine

import cats.data.NonEmptyList
import cromwell.CromwellTestKitSpec
import cromwell.core.{WorkflowId, WorkflowSourceFilesCollection}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore._
import cromwell.services.metadata.MetadataQuery
import cromwell.services.metadata.MetadataService.{GetMetadataQueryAction, MetadataLookupResponse}
import cromwell.services.metadata.impl.ReadMetadataActor
import cromwell.util.EncryptionSpec
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.Matchers

import scala.concurrent.duration._
import scala.language.postfixOps

class WorkflowStoreActorSpec extends CromwellTestKitSpec with Matchers {
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

  private def prettyOptions(workflowSourceFiles: WorkflowSourceFilesCollection): WorkflowSourceFilesCollection = {
    import spray.json._
    workflowSourceFiles.copyOptions(workflowSourceFiles.workflowOptionsJson.parseJson.prettyPrint)
  }

  "The WorkflowStoreActor" should {
    "return an ID for a submitted workflow" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store, CromwellTestKitSpec.ServiceRegistryActorInstance))
      storeActor ! SubmitWorkflow(helloWorldSourceFiles)
      expectMsgType[WorkflowSubmittedToStore](10 seconds)
    }

    "return 3 IDs for a batch submission of 3" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store, CromwellTestKitSpec.ServiceRegistryActorInstance))
      storeActor ! BatchSubmitWorkflows(NonEmptyList.of(helloWorldSourceFiles, helloWorldSourceFiles, helloWorldSourceFiles))
      expectMsgPF(10 seconds) {
        case WorkflowsBatchSubmittedToStore(ids) => ids.toList.size shouldBe 3
      }
    }

    "fetch exactly N workflows" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store, CromwellTestKitSpec.ServiceRegistryActorInstance))
      storeActor ! BatchSubmitWorkflows(NonEmptyList.of(helloWorldSourceFiles, helloWorldSourceFiles, helloWorldSourceFiles))
      val insertedIds = expectMsgType[WorkflowsBatchSubmittedToStore](10 seconds).workflowIds.toList


      storeActor ! FetchRunnableWorkflows(2)
      expectMsgPF(10 seconds) {
        case NewWorkflowsToStart(workflowNel) =>
          workflowNel.toList.size shouldBe 2
          checkDistinctIds(workflowNel.toList) shouldBe true
          workflowNel map {
            case WorkflowToStart(id, sources, state) =>
              insertedIds.contains(id) shouldBe true
              sources shouldBe prettyOptions(helloWorldSourceFiles)
              state shouldBe WorkflowStoreState.Submitted
          }
      }
    }

    "fetch encrypted and cleared workflow options" in {
      EncryptionSpec.assumeAes256Cbc()

      val optionedSourceFiles = HelloWorld.asWorkflowSources(workflowOptions =
        s"""|{
            |  "key": "value",
            |  "refresh_token": "it's a secret"
            |}
            |""".stripMargin)


      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store, CromwellTestKitSpec.ServiceRegistryActorInstance))
      val readMetadataActor = system.actorOf(ReadMetadataActor.props())
      storeActor ! BatchSubmitWorkflows(NonEmptyList.of(optionedSourceFiles))
      val insertedIds = expectMsgType[WorkflowsBatchSubmittedToStore](10 seconds).workflowIds.toList

      storeActor ! FetchRunnableWorkflows(1)
      expectMsgPF(10 seconds) {
        case NewWorkflowsToStart(workflowNel) =>
          workflowNel.toList.size should be(1)
          checkDistinctIds(workflowNel.toList) should be(true)
          workflowNel.toList.foreach {
            case WorkflowToStart(id, sources, state) =>
              insertedIds.contains(id) should be(true)
              sources.wdlSource should be(optionedSourceFiles.wdlSource)
              sources.inputsJson should be(optionedSourceFiles.inputsJson)
              state should be(WorkflowStoreState.Submitted)

              import spray.json._

              val encryptedJsObject = sources.workflowOptionsJson.parseJson.asJsObject
              encryptedJsObject.fields.keys should contain theSameElementsAs Seq("key", "refresh_token")
              encryptedJsObject.fields("key") should be(JsString("value"))
              encryptedJsObject.fields("refresh_token").asJsObject.fields.keys should contain theSameElementsAs
                Seq("iv", "ciphertext")

              readMetadataActor ! GetMetadataQueryAction(MetadataQuery.forWorkflow(id))
              expectMsgPF(10 seconds) {
                case MetadataLookupResponse(_, eventList) =>
                  val optionsEvent = eventList.find(_.key.key == "submittedFiles:options").get
                  val clearedJsObject = optionsEvent.value.get.value.parseJson.asJsObject
                  clearedJsObject.fields.keys should contain theSameElementsAs Seq("key", "refresh_token")
                  clearedJsObject.fields("key") should be(JsString("value"))
                  clearedJsObject.fields("refresh_token") should be(JsString("cleared"))
              }
          }
      }
    }

    "return only the remaining workflows if N is larger than size" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store, CromwellTestKitSpec.ServiceRegistryActorInstance))
      storeActor ! BatchSubmitWorkflows(NonEmptyList.of(helloWorldSourceFiles, helloWorldSourceFiles, helloWorldSourceFiles))
      val insertedIds = expectMsgType[WorkflowsBatchSubmittedToStore](10 seconds).workflowIds.toList


      storeActor ! FetchRunnableWorkflows(100)
      expectMsgPF(10 seconds) {
        case NewWorkflowsToStart(workflowNel) =>
          workflowNel.toList.size shouldBe 3
          checkDistinctIds(workflowNel.toList) shouldBe true
          workflowNel map {
            case WorkflowToStart(id, sources, state) =>
              insertedIds.contains(id) shouldBe true
              sources shouldBe prettyOptions(helloWorldSourceFiles)
              state shouldBe WorkflowStoreState.Submitted
          }
      }
    }

    "remove workflows which exist" in {
      val store = new InMemoryWorkflowStore
      val storeActor = system.actorOf(WorkflowStoreActor.props(store, CromwellTestKitSpec.ServiceRegistryActorInstance))
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
      val storeActor = system.actorOf(WorkflowStoreActor.props(store, CromwellTestKitSpec.ServiceRegistryActorInstance))
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
