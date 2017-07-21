package cromwell.jobstore

import cromwell.CromwellTestKitWordSpec
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.{JobOutput, WorkflowId}
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore.JobStoreServiceSpec._
import cromwell.services.SingletonServicesStore
import org.scalatest.Matchers
import org.specs2.mock.Mockito
import wdl4s.parser.WdlParser.Ast
import wdl4s.wdl.types.WdlStringType
import wdl4s.wdl.values.WdlString
import wdl4s.wdl._

import scala.concurrent.duration._
import scala.language.postfixOps

object JobStoreServiceSpec {
  val MaxWait = 5 seconds
  val EmptyExpression = WdlExpression.fromString(""" "" """)
}

class JobStoreServiceSpec extends CromwellTestKitWordSpec with Matchers with Mockito {

  "JobStoreService" should {
    "work" in {
      lazy val jobStore: JobStore = new SqlJobStore(SingletonServicesStore.databaseInterface)
      val jobStoreService = system.actorOf(JobStoreActor.props(jobStore))

      val workflowId = WorkflowId.randomId()
      val successCall = mock[WdlTaskCall]
      successCall.fullyQualifiedName returns "foo.bar"
      val mockTask = mock[WdlTask]
      mockTask.outputs returns Seq(TaskOutput("baz", WdlStringType, EmptyExpression, mock[Ast], Option(mockTask)))
      successCall.task returns mockTask

      val successKey = BackendJobDescriptorKey(successCall, None, 1).toJobStoreKey(workflowId)

      jobStoreService ! QueryJobCompletion(successKey, mockTask.outputs)
      expectMsgType[JobNotComplete.type](MaxWait)

      val outputs = Map("baz" -> JobOutput(WdlString("qux")))

      jobStoreService ! RegisterJobCompleted(successKey, JobResultSuccess(Option(0), outputs))
      expectMsgType[JobStoreWriteSuccess](MaxWait)

      jobStoreService ! QueryJobCompletion(successKey, mockTask.outputs)
      expectMsgPF(MaxWait) {
        case JobComplete(JobResultSuccess(Some(0), os)) if os == outputs =>
      }

      val failureCall = mock[WdlTaskCall]
      failureCall.fullyQualifiedName returns "baz.qux"
      val failureKey = BackendJobDescriptorKey(failureCall, None, 1).toJobStoreKey(workflowId)

      jobStoreService ! QueryJobCompletion(failureKey, mockTask.outputs)
      expectMsgType[JobNotComplete.type](MaxWait)

      jobStoreService ! RegisterJobCompleted(failureKey, JobResultFailure(Option(11), new IllegalArgumentException("Insufficient funds"), retryable = false))
      expectMsgType[JobStoreWriteSuccess](MaxWait)

      jobStoreService ! QueryJobCompletion(failureKey, mockTask.outputs)
      expectMsgPF(MaxWait) {
        case JobComplete(JobResultFailure(Some(11), _, false)) =>
      }

      jobStoreService ! RegisterWorkflowCompleted(workflowId)
      expectMsgType[JobStoreWriteSuccess](MaxWait)
    }
  }
}
