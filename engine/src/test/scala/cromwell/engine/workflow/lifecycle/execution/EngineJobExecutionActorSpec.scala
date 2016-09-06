package cromwell.engine.workflow.lifecycle.execution

import java.util.UUID

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.backend.BackendJobExecutionActor._
import cromwell.backend._
import cromwell.core.callcaching.CallCachingOff
import cromwell.core.{ExecutionStore, JobOutputs, OutputStore, WorkflowId}
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor.BackendJobPreparationFailed
import cromwell.jobstore.JobStoreActor.{JobComplete, JobNotComplete}
import cromwell.jobstore.{JobResultFailure, JobResultSuccess, JobStoreActor, SqlJobStore, Pending => _}
import cromwell.services.SingletonServicesStore
import cromwell.util.SampleWdl
import cromwell.{CromwellTestkitSpec, EmptyCallCacheReadActor}
import org.scalatest.{BeforeAndAfterAll, Matchers}
import org.specs2.mock.Mockito
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions}
import wdl4s.{Call, RuntimeAttributes, Task}

import scala.concurrent.duration._
import scala.language.postfixOps

class EngineJobExecutionActorSpec extends CromwellTestkitSpec with Matchers with WorkflowDescriptorBuilder with Mockito with BeforeAndAfterAll {

  override implicit val actorSystem = system

  private val workflowId: WorkflowId = WorkflowId.randomId()
  val descriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.HelloWorld.asWorkflowSources())
  val messageProbe = TestProbe()
  val replyToProbe = TestProbe()
  val parentProbe = TestProbe()
  val mockBackendProps = Props(new Actor {
    def receive = {
      case x => messageProbe.ref forward x
    }
  })
  val factory = new BackendLifecycleActorFactory {
    override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                  calls: Seq[Call],
                                                  serviceRegistryActor: ActorRef): Option[Props] = None
    override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                        initializationData: Option[BackendInitializationData],
                                        serviceRegistryActor: ActorRef): Props = mockBackendProps

    override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor, jobKey: BackendJobDescriptorKey, initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {
      NoFunctions
    }
  }

  def buildEJEA(restarting: Boolean) = {

    val task = mock[Task]
    task.declarations returns Seq.empty
    task.runtimeAttributes returns RuntimeAttributes(Map.empty)

    val jobStore = new SqlJobStore(SingletonServicesStore.databaseInterface)

    val callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props)

    new TestFSMRef[EngineJobExecutionActorState, EJEAData, EngineJobExecutionActor](system, EngineJobExecutionActor.props(
      replyTo = replyToProbe.ref,
      BackendJobDescriptorKey(Call(None, "foo.bar", task, Set.empty, Map.empty, None), None, 1),
      WorkflowExecutionActorData(descriptor, ExecutionStore(Map.empty), Map.empty, OutputStore(Map.empty)),
      factory,
      None,
      restarting = restarting,
      serviceRegistryActor = CromwellTestkitSpec.ServiceRegistryActorInstance,
      jobStoreActor = system.actorOf(JobStoreActor.props(jobStore)),
      callCacheReadActor = callCacheReadActor,
      backendName = "NOT USED",
      callCachingMode = CallCachingOff
    ), parentProbe.ref, s"EngineJobExecutionActorSpec-$workflowId")
  }

  "EngineJobExecutionActorSpec" should {
    "send a Job SucceededResponse if the job is already complete and successful" in {
      val ejea = buildEJEA(restarting = true)
      ejea.setState(CheckingJobStore)
      val returnCode: Option[Int] = Option(0)
      val jobOutputs: JobOutputs = Map.empty

      ejea ! JobComplete(JobResultSuccess(returnCode, jobOutputs))

      replyToProbe.expectMsgPF(awaitTimeout) {
        case response: SucceededResponse =>
          response.returnCode shouldBe returnCode
          response.jobOutputs shouldBe jobOutputs
      }

      ejea.stop()
    }

    "send a Job FailedNonRetryableResponse if the job is already complete and failed" in {
      val ejea = buildEJEA(restarting = true)
      ejea.setState(CheckingJobStore)
      val returnCode: Option[Int] = Option(1)
      val reason: Throwable = new Exception("something horrible happened...")

      ejea ! JobComplete(JobResultFailure(returnCode, reason, retryable = false))

      replyToProbe.expectMsgPF(awaitTimeout) {
        case response: FailedNonRetryableResponse =>
          response.returnCode shouldBe returnCode
          response.throwable shouldBe reason
      }

      ejea.stop()
    }

    "send a RecoverJobCommand to the backend if the job is not complete and EJEA is in restarting mode" in {
      val ejea = buildEJEA(restarting = true)
      ejea.setState(CheckingJobStore)
      ejea ! JobNotComplete

      messageProbe.expectMsg(10 seconds, "expecting RecoverJobCommand", RecoverJobCommand)
      replyToProbe.expectMsgType[JobRunning]

      ejea.stop()
    }

    // This test is bit redundant with the one above, however this one tests the behavior when getting a Restart message
    // whereas the test above tests the behavior when getting a JobNotComplete response from JobReader
    "send a RecoverJobCommand to the backend when getting an Execute message if the Job is not complete and EJEA is in restarting mode" in {
      val ejea = buildEJEA(restarting = true)
      ejea.setState(Pending)

      ejea ! EngineJobExecutionActor.Execute

      replyToProbe.expectMsgType[JobRunning](awaitTimeout)
      messageProbe.expectMsg(awaitTimeout, "RecoverJobCommand", RecoverJobCommand)

      ejea.stop()
    }

    "send an ExecuteJobCommand to the backend when getting an Execute message if the Job is not complete and EJEA is NOT in restarting mode" in {
      val ejea = buildEJEA(restarting = false)
      ejea.setState(Pending)

      ejea ! EngineJobExecutionActor.Execute

      messageProbe.expectMsg(awaitTimeout, ExecuteJobCommand)
      replyToProbe.expectMsgType[JobRunning](awaitTimeout)

      ejea.stop()
    }

    "forward JobPreparationFailed to its parent and die" in {
      val deathWatch = TestProbe()
      val jobKey = mock[BackendJobDescriptorKey]

      val ejea = buildEJEA(restarting = true)
      ejea.setState(PreparingJob)

      val failed = BackendJobPreparationFailed(jobKey, new scala.Exception("Failed"))
      deathWatch watch ejea

      ejea ! failed
      replyToProbe.expectMsg(awaitTimeout, failed)

      deathWatch.expectTerminated(ejea)
    }

    "forward responds from backend to its parent" in {
      val deathWatch = TestProbe()

      def jobKey = {
        val jk = mock[BackendJobDescriptorKey]
        val scope = mock[Call]
        jk.scope returns scope
        scope.fullyQualifiedName returns "blah." + UUID.randomUUID().toString // Make sure job keys are unique.
        jk.index returns None
        jk
      }

      def ejeaInRunningState() = {
        val ejea = buildEJEA(restarting = true)
        ejea.setState(stateName = RunningJob, stateData = ResponsePendingData(null, Props.empty))
        ejea
      }

      val successResponse = SucceededResponse(jobKey, Option(0), Map.empty)
      val failureRetryableResponse = FailedRetryableResponse(jobKey, new Exception("Failed"), Option(0))
      val failureNonRetryableResponse = FailedNonRetryableResponse(jobKey, new Exception("Failed"), Option(0))
      val abortedResponse = FailedNonRetryableResponse(jobKey, new Exception("Failed"), Option(0))

      val ejea1 = ejeaInRunningState()
      deathWatch watch ejea1
      ejea1 ! successResponse
      replyToProbe.expectMsg(awaitTimeout, successResponse)
      deathWatch.expectTerminated(ejea1)

      val ejea2 = ejeaInRunningState()
      deathWatch watch ejea2
      ejea2 ! failureRetryableResponse
      replyToProbe.expectMsg(awaitTimeout, failureRetryableResponse)
      deathWatch.expectTerminated(ejea2)

      val ejea3 = ejeaInRunningState()
      deathWatch watch ejea3
      ejea3 ! failureNonRetryableResponse
      replyToProbe.expectMsg(awaitTimeout, failureNonRetryableResponse)
      deathWatch.expectTerminated(ejea3)

      val ejea4 = ejeaInRunningState()
      deathWatch watch ejea4
      ejea4 ! abortedResponse
      replyToProbe.expectMsg(awaitTimeout, abortedResponse)
      deathWatch.expectTerminated(ejea4)
    }
  }

  override def afterAll(): Unit = {
    system.terminate()
  }
}
