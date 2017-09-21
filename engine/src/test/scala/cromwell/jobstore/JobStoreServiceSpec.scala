package cromwell.jobstore

import cromwell.CromwellTestKitWordSpec
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.{JobOutput, WorkflowId}
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore.JobStoreServiceSpec._
import cromwell.services.EngineServicesStore
import cromwell.util.WomMocks
import org.scalatest.Matchers
import org.specs2.mock.Mockito
import wdl4s.wdl.types.WdlStringType
import wdl4s.wdl.values.WdlString
import wdl4s.wom.callable.Callable.OutputDefinition
import wdl4s.wom.expression.PlaceholderWomExpression

import scala.concurrent.duration._
import scala.language.postfixOps

object JobStoreServiceSpec {
  val MaxWait = 5 seconds
  val EmptyExpression = PlaceholderWomExpression(Set.empty, WdlStringType)
}

class JobStoreServiceSpec extends CromwellTestKitWordSpec with Matchers with Mockito {

  "JobStoreService" should {
    "work" in {
      lazy val jobStore: JobStore = new SqlJobStore(EngineServicesStore.engineDatabaseInterface)
      val jobStoreService = system.actorOf(JobStoreActor.props(jobStore))

      val workflowId = WorkflowId.randomId()
      val mockTask = WomMocks.mockTaskDefinition("bar")
      .copy(outputs = Set(OutputDefinition("baz", WdlStringType, EmptyExpression))) 
      val successCall = WomMocks.mockTaskCall("bar", definition = mockTask)

      val successKey = BackendJobDescriptorKey(successCall, None, 1).toJobStoreKey(workflowId)

      jobStoreService ! QueryJobCompletion(successKey, mockTask.outputs.toSeq)
      expectMsgType[JobNotComplete.type](MaxWait)

      val outputs = Map("baz" -> JobOutput(WdlString("qux")))

      jobStoreService ! RegisterJobCompleted(successKey, JobResultSuccess(Option(0), outputs))
      expectMsgType[JobStoreWriteSuccess](MaxWait)

      jobStoreService ! QueryJobCompletion(successKey, mockTask.outputs.toSeq)
      expectMsgPF(MaxWait) {
        case JobComplete(JobResultSuccess(Some(0), os)) if os == outputs =>
      }

      val failureCall = WomMocks.mockTaskCall("qux")
      val failureKey = BackendJobDescriptorKey(failureCall, None, 1).toJobStoreKey(workflowId)

      jobStoreService ! QueryJobCompletion(failureKey, mockTask.outputs.toSeq)
      expectMsgType[JobNotComplete.type](MaxWait)

      jobStoreService ! RegisterJobCompleted(failureKey, JobResultFailure(Option(11), new IllegalArgumentException("Insufficient funds"), retryable = false))
      expectMsgType[JobStoreWriteSuccess](MaxWait)

      jobStoreService ! QueryJobCompletion(failureKey, mockTask.outputs.toSeq)
      expectMsgPF(MaxWait) {
        case JobComplete(JobResultFailure(Some(11), _, false)) =>
      }

      jobStoreService ! RegisterWorkflowCompleted(workflowId)
      expectMsgType[JobStoreWriteSuccess](MaxWait)
    }
  }
}
