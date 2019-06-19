package cromwell.subworkflowstore

import akka.testkit.TestProbe
import cromwell.core.ExecutionIndex._
import cromwell.core.{JobKey, WorkflowId, WorkflowSourceFilesWithoutImports}
import cromwell.database.sql.tables.SubWorkflowStoreEntry
import cromwell.engine.MockCromwellTerminator
import cromwell.engine.workflow.CoordinatedWorkflowStoreBuilder
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.SubmitWorkflow
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.WorkflowSubmittedToStore
import cromwell.engine.workflow.workflowstore._
import cromwell.services.EngineServicesStore
import cromwell.subworkflowstore.SubWorkflowStoreActor._
import cromwell.subworkflowstore.SubWorkflowStoreSpec._
import cromwell.util.WomMocks
import cromwell.{CromwellTestKitSpec, CromwellTestKitWordSpec}
import mouse.all._
import org.scalatest.Matchers
import org.specs2.mock.Mockito
import wdl.draft2.model.WdlExpression
import wom.graph.{GraphNode, WomIdentifier}

import scala.concurrent.duration._
import scala.language.postfixOps

object SubWorkflowStoreSpec {
  val MaxWait = 5 seconds
  val EmptyExpression = WdlExpression.fromString(""" "" """)
}

class SubWorkflowStoreSpec extends CromwellTestKitWordSpec with CoordinatedWorkflowStoreBuilder with Matchers with Mockito {
  "SubWorkflowStore" should {
    "work" in {
      lazy val subWorkflowStore = new SqlSubWorkflowStore(EngineServicesStore.engineDatabaseInterface)
      val subWorkflowStoreService = system.actorOf(SubWorkflowStoreActor.props(subWorkflowStore))

      lazy val workflowStore = SqlWorkflowStore(EngineServicesStore.engineDatabaseInterface)
      val workflowHeartbeatConfig = WorkflowHeartbeatConfig(CromwellTestKitSpec.DefaultConfig)
      val workflowStoreService = system.actorOf(
        WorkflowStoreActor.props(
          workflowStore,
          workflowStore |> access,
          TestProbe("ServiceRegistryProbe-Work").ref,
          MockCromwellTerminator,
          abortAllJobsOnTerminate = false,
          workflowHeartbeatConfig
        ),
        "WorkflowStoreActor-Work"
      )

      val parentWorkflowId = WorkflowId.randomId()
      val subWorkflowId = WorkflowId.randomId()
      val subSubWorkflowId = WorkflowId.randomId()
      val call = WomMocks.mockTaskCall(WomIdentifier("bar", "foo.bar"))
      val jobKey = new JobKey {
        override def node: GraphNode = call
        override def index: Option[Int] = None
        override def attempt: Int = 0
        override def tag: String = "foobar"
      }

      workflowStoreService ! SubmitWorkflow(WorkflowSourceFilesWithoutImports(
        workflowSource = Option(""),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptionsJson = "{}",
        labelsJson = "{}",
        warnings = Vector.empty)
      )
      val rootWorkflowId = expectMsgType[WorkflowSubmittedToStore](10 seconds).workflowId

      // Query for non existing sub workflow
      subWorkflowStoreService ! QuerySubWorkflow(parentWorkflowId, jobKey)
      expectMsgType[SubWorkflowNotFound](MaxWait)
      
      // Register sub workflow
      subWorkflowStoreService ! RegisterSubWorkflow(rootWorkflowId, parentWorkflowId, jobKey, subWorkflowId)
      expectMsgType[SubWorkflowStoreRegisterSuccess](MaxWait)

      // Query for sub workflow
      subWorkflowStoreService ! QuerySubWorkflow(parentWorkflowId, jobKey)
      val subWorkflowEntry = SubWorkflowStoreEntry(Option(0), parentWorkflowId.toString, jobKey.node.fullyQualifiedName, jobKey.index.fromIndex, jobKey.attempt, subWorkflowId.toString, Some(0))
      expectMsg[SubWorkflowFound](SubWorkflowFound(subWorkflowEntry))

      // Register sub sub workflow
      subWorkflowStoreService ! RegisterSubWorkflow(rootWorkflowId, subWorkflowId, jobKey, subSubWorkflowId)
      expectMsgType[SubWorkflowStoreRegisterSuccess](MaxWait)

      // Query for sub sub workflow
      subWorkflowStoreService ! QuerySubWorkflow(subWorkflowId, jobKey)
      val subSubWorkflowEntry = SubWorkflowStoreEntry(Option(0), subWorkflowId.toString, jobKey.node.fullyQualifiedName, jobKey.index.fromIndex, jobKey.attempt, subSubWorkflowId.toString, Some(1))
      expectMsg[SubWorkflowFound](SubWorkflowFound(subSubWorkflowEntry))
      
      // Delete root workflow
      subWorkflowStoreService ! WorkflowComplete(rootWorkflowId)

      // Verify that everything is gone (eventually!)
      eventually {
        subWorkflowStoreService ! QuerySubWorkflow(rootWorkflowId, jobKey)
        expectMsgType[SubWorkflowNotFound](MaxWait)
      }

      eventually {
        subWorkflowStoreService ! QuerySubWorkflow(parentWorkflowId, jobKey)
        expectMsgType[SubWorkflowNotFound](MaxWait)
      }

      eventually {
        subWorkflowStoreService ! QuerySubWorkflow(subWorkflowId, jobKey)
        expectMsgType[SubWorkflowNotFound](MaxWait)
      }
    }
  }
}
