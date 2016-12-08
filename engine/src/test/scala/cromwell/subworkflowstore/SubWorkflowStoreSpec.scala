package cromwell.subworkflowstore

import cromwell.CromwellTestKitSpec
import cromwell.core.{JobKey, WorkflowId, WorkflowSourceFilesWithoutImports}
import cromwell.services.SingletonServicesStore
import cromwell.subworkflowstore.SubWorkflowStoreActor._
import org.scalatest.Matchers
import org.specs2.mock.Mockito
import wdl4s.{TaskCall, WdlExpression}
import cromwell.core.ExecutionIndex._

import scala.concurrent.duration._
import SubWorkflowStoreSpec._
import akka.testkit.TestProbe
import cromwell.database.sql.tables.SubWorkflowStoreEntry
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.{SubmitWorkflow, WorkflowSubmittedToStore}
import cromwell.engine.workflow.workflowstore.{SqlWorkflowStore, WorkflowStoreActor}

import scala.language.postfixOps

object SubWorkflowStoreSpec {
  val MaxWait = 5 seconds
  val EmptyExpression = WdlExpression.fromString(""" "" """)
}

class SubWorkflowStoreSpec extends CromwellTestKitSpec with Matchers with Mockito {
  "SubWorkflowStore" should {
    "work" in {
      lazy val subWorkflowStore = new SqlSubWorkflowStore(SingletonServicesStore.databaseInterface)
      val subWorkflowStoreService = system.actorOf(SubWorkflowStoreActor.props(subWorkflowStore))

      lazy val workflowStore = SqlWorkflowStore(SingletonServicesStore.databaseInterface)
      val workflowStoreService = system.actorOf(WorkflowStoreActor.props(workflowStore, TestProbe().ref))

      val parentWorkflowId = WorkflowId.randomId()
      val subWorkflowId = WorkflowId.randomId()
      val subSubWorkflowId = WorkflowId.randomId()
      val call = mock[TaskCall]
      call.fullyQualifiedName returns "foo.bar"
      val jobKey = new JobKey {
        override def scope = call
        override def index: Option[Int] = None
        override def attempt: Int = 0
        override def tag: String = "foobar"
      }

      workflowStoreService ! SubmitWorkflow(WorkflowSourceFilesWithoutImports("", "{}", "{}"))
      val rootWorkflowId = expectMsgType[WorkflowSubmittedToStore](10 seconds).workflowId

      // Query for non existing sub workflow
      subWorkflowStoreService ! QuerySubWorkflow(parentWorkflowId, jobKey)
      expectMsgType[SubWorkflowNotFound](MaxWait)
      
      // Register sub workflow
      subWorkflowStoreService ! RegisterSubWorkflow(rootWorkflowId, parentWorkflowId, jobKey, subWorkflowId)
      expectMsgType[SubWorkflowStoreRegisterSuccess](MaxWait)

      // Query for sub workflow
      subWorkflowStoreService ! QuerySubWorkflow(parentWorkflowId, jobKey)
      val subWorkflowEntry = SubWorkflowStoreEntry(Option(0), parentWorkflowId.toString, jobKey.scope.fullyQualifiedName, jobKey.index.fromIndex, jobKey.attempt, subWorkflowId.toString, Some(0))
      expectMsg[SubWorkflowFound](SubWorkflowFound(subWorkflowEntry))

      // Register sub sub workflow
      subWorkflowStoreService ! RegisterSubWorkflow(rootWorkflowId, subWorkflowId, jobKey, subSubWorkflowId)
      expectMsgType[SubWorkflowStoreRegisterSuccess](MaxWait)

      // Query for sub sub workflow
      subWorkflowStoreService ! QuerySubWorkflow(subWorkflowId, jobKey)
      val subSubWorkflowEntry = SubWorkflowStoreEntry(Option(0), subWorkflowId.toString, jobKey.scope.fullyQualifiedName, jobKey.index.fromIndex, jobKey.attempt, subSubWorkflowId.toString, Some(1))
      expectMsg[SubWorkflowFound](SubWorkflowFound(subSubWorkflowEntry))
      
      // Delete root workflow
      subWorkflowStoreService ! WorkflowComplete(rootWorkflowId)
      expectMsgType[SubWorkflowStoreCompleteSuccess](MaxWait)

      // Verify that everything is gone
      subWorkflowStoreService ! QuerySubWorkflow(rootWorkflowId, jobKey)
      expectMsgType[SubWorkflowNotFound](MaxWait)

      subWorkflowStoreService ! QuerySubWorkflow(parentWorkflowId, jobKey)
      expectMsgType[SubWorkflowNotFound](MaxWait)

      subWorkflowStoreService ! QuerySubWorkflow(subWorkflowId, jobKey)
      expectMsgType[SubWorkflowNotFound](MaxWait)
    }
  }
}
