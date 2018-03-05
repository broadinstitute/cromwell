package cromwell.jobstore

import cromwell.CromwellTestKitWordSpec
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.WorkflowId
import cromwell.jobstore.JobStoreActor._
import cromwell.jobstore.JobStoreServiceSpec._
import cromwell.services.EngineServicesStore
import cromwell.util.WomMocks
import org.scalatest.Matchers
import org.specs2.mock.Mockito
import wom.callable.Callable.OutputDefinition
import wom.expression.PlaceholderWomExpression
import wom.graph.WomIdentifier
import wom.types.WomStringType
import wom.values.WomString

import scala.concurrent.duration._
import scala.language.postfixOps

object JobStoreServiceSpec {
  val MaxWait = 5 seconds
  val EmptyExpression = PlaceholderWomExpression(Set.empty, WomStringType)
}

class JobStoreServiceSpec extends CromwellTestKitWordSpec with Matchers with Mockito {

  "JobStoreService" should {
    "work" in {
      lazy val jobStore: JobStore = new SqlJobStore(EngineServicesStore.engineDatabaseInterface)
      val jobStoreService = system.actorOf(JobStoreActor.props(jobStore, dummyServiceRegistryActor))

      val workflowId = WorkflowId.randomId()
      val mockTask = WomMocks.mockTaskDefinition("bar")
      .copy(outputs = List(OutputDefinition("baz", WomStringType, EmptyExpression)))
      val successCall = WomMocks.mockTaskCall(WomIdentifier("bar"), definition = mockTask)

      val successKey = BackendJobDescriptorKey(successCall, None, 1).toJobStoreKey(workflowId)

      jobStoreService ! QueryJobCompletion(successKey, mockTask.outputs map WomMocks.mockOutputPort)
      expectMsgType[JobNotComplete.type](MaxWait)

      val outputs = WomMocks.mockOutputExpectations(Map("baz" -> WomString("qux")))

      jobStoreService ! RegisterJobCompleted(successKey, JobResultSuccess(Option(0), outputs))
      expectMsgType[JobStoreWriteSuccess](MaxWait)

      jobStoreService ! QueryJobCompletion(successKey, mockTask.outputs map WomMocks.mockOutputPort)
      expectMsgPF(MaxWait) {
        case JobComplete(JobResultSuccess(Some(0), os)) if os == outputs =>
      }

      val failureCall = WomMocks.mockTaskCall(WomIdentifier("qux"))
      val failureKey = BackendJobDescriptorKey(failureCall, None, 1).toJobStoreKey(workflowId)

      jobStoreService ! QueryJobCompletion(failureKey, mockTask.outputs map WomMocks.mockOutputPort)
      expectMsgType[JobNotComplete.type](MaxWait)

      jobStoreService ! RegisterJobCompleted(failureKey, JobResultFailure(Option(11), new IllegalArgumentException("Insufficient funds"), retryable = false))
      expectMsgType[JobStoreWriteSuccess](MaxWait)

      jobStoreService ! QueryJobCompletion(failureKey, mockTask.outputs map WomMocks.mockOutputPort)
      expectMsgPF(MaxWait) {
        case JobComplete(JobResultFailure(Some(11), _, false)) =>
      }

      jobStoreService ! RegisterWorkflowCompleted(workflowId)
      expectMsgType[JobStoreWriteSuccess](MaxWait)
    }
  }
}
