package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.TestProbe
import cromwell.backend._
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData
import org.mockito.Mockito._
import org.specs2.mock.Mockito
import wdl4s.Declaration
import wdl4s.values.WdlValue

import scala.concurrent.duration.FiniteDuration
import scala.util.Try

class JobPreparationTestHelper(implicit val system: ActorSystem) extends Mockito {
  val executionData = mock[WorkflowExecutionActorData]
  val workflowDescriptor = mock[EngineWorkflowDescriptor]
  val backendDescriptor = mock[BackendWorkflowDescriptor]
  workflowDescriptor.backendDescriptor returns backendDescriptor
  executionData.workflowDescriptor returns workflowDescriptor
  val jobKey = mock[BackendJobDescriptorKey]
  val serviceRegistryProbe = TestProbe()
  val ioActor = TestProbe()

  def buildJobPreparationMock(
                               backpressureTimeoutValue: FiniteDuration,
                               noResponsTimeoutValue: FiniteDuration,
                               dockerHashCredentials: List[Any],
                               dockerHashingActor: ActorRef,
                               inputsAndAttributes: Try[(Map[Declaration, WdlValue], Map[wdl4s.LocallyQualifiedName, WdlValue])]
                             ) = {
    val factory = mock[BackendLifecycleActorFactory]
    when(factory.dockerHashCredentials(None)).thenReturn(dockerHashCredentials)
    when(factory.jobExecutionActorProps(
      any[BackendJobDescriptor],
      any[Option[BackendInitializationData]],
      any[ActorRef],
      any[ActorRef],
      any[Option[ActorRef]]
    )).thenReturn(Props.empty)

    Props(new JobPreparationActor(
      executionData,
      jobKey,
      factory,
      dockerHashingActor,
      None,
      ioActor.ref,
      serviceRegistryProbe.ref,
      None
    ) {
      override protected def backpressureTimeout: FiniteDuration = backpressureTimeoutValue
      override def evaluateInputsAndAttributes = inputsAndAttributes
      override protected lazy val noResponsTimeout = noResponsTimeoutValue
    })
  }
}
