package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{Actor, Props}
import akka.testkit.{EventFilter, TestActorRef, TestDuration, TestProbe}
import com.typesafe.config.ConfigFactory
import cromwell._
import cromwell.backend.AllBackendInitializationData
import cromwell.core.{SimpleIoActor, WorkflowId}
import cromwell.engine.backend.{BackendConfigurationEntry, BackendSingletonCollection, CromwellBackends}
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{ExecuteWorkflowCommand, WorkflowExecutionFailedResponse}
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor
import cromwell.services.ServiceRegistryActor
import cromwell.services.metadata.MetadataService
import cromwell.util.SampleWdl
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{Await, Promise}

class WorkflowExecutionActorSpec extends CromwellTestKitSpec with FlatSpecLike with Matchers with BeforeAndAfter with WorkflowDescriptorBuilder {

  override implicit val actorSystem = system
  implicit val DefaultDuration = 20.seconds.dilated

  def mockServiceRegistryActor = TestActorRef(new Actor {
    override def receive = {
      case _ => // No action
    }
  })

  val MockBackendName = "Mock"
  val MockBackendSingletonCollection = BackendSingletonCollection(Map(MockBackendName -> None))

  val stubbedConfig = ConfigFactory.load().getConfig("backend.providers.Mock").getConfig("config")

  val runtimeSection =
    """
      |runtime {
      | backend: "Mock"
      |}
    """.stripMargin

  behavior of "WorkflowExecutionActor"

  it should "allow a backend to tell it to retry... up to a point" in {
    import MetadataWatchActor.metadataKeyAttemptChecker
    val metadataSuccessPromise = Promise[Unit]()
    val requiredMetadataMatchers: Seq[MetadataWatchActor.Matcher] = List(
      MetadataWatchActor.JobKeyMetadataKeyAndValueContainStringMatcher(metadataKeyAttemptChecker(1), "executionStatus", "RetryableFailure"),
      MetadataWatchActor.JobKeyMetadataKeyAndValueContainStringMatcher(metadataKeyAttemptChecker(2), "executionStatus", "RetryableFailure"),
      MetadataWatchActor.JobKeyMetadataKeyAndValueContainStringMatcher(metadataKeyAttemptChecker(3), "executionStatus", "Failed")
    )
    val metadataWatcherActor = TestActorRef[MetadataWatchActor](Props(MetadataWatchActor(metadataSuccessPromise, requiredMetadataMatchers: _*)))
    val serviceRegistryActor = system.actorOf(ServiceRegistryActor.props(ConfigFactory.load(), overrides = Map(MetadataService.MetadataServiceName -> metadataWatcherActor.props)))
    val jobStoreActor = system.actorOf(AlwaysHappyJobStoreActor.props)
    val ioActor = system.actorOf(SimpleIoActor.props)
    val subWorkflowStoreActor = system.actorOf(AlwaysHappySubWorkflowStoreActor.props)
    val jobTokenDispenserActor = system.actorOf(JobExecutionTokenDispenserActor.props)
    val MockBackendConfigEntry = BackendConfigurationEntry(
      name = "Mock",
      lifecycleActorFactoryClass = "cromwell.engine.backend.mock.RetryableBackendLifecycleActorFactory",
      stubbedConfig
    )
    CromwellBackends.initBackends(List(MockBackendConfigEntry))

    val workflowId = WorkflowId.randomId()
    val engineWorkflowDescriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.HelloWorld.asWorkflowSources(runtime = runtimeSection))
    val callCacheReadActor = TestProbe()
    val callCacheWriteActor = TestProbe()
    val dockerHashActor = TestProbe()

    val weaSupervisor = TestProbe()
    val workflowExecutionActor = TestActorRef(
      props = WorkflowExecutionActor.props(engineWorkflowDescriptor, ioActor, serviceRegistryActor, jobStoreActor, subWorkflowStoreActor,
        callCacheReadActor.ref, callCacheWriteActor.ref, dockerHashActor.ref, jobTokenDispenserActor, MockBackendSingletonCollection, AllBackendInitializationData.empty, restarting = false),
      name = "WorkflowExecutionActor",
      supervisor = weaSupervisor.ref)

    EventFilter.info(pattern = "Starting calls: wf_hello.hello", occurrences = 3).intercept {
      workflowExecutionActor ! ExecuteWorkflowCommand
    }

    weaSupervisor.expectMsgClass(classOf[WorkflowExecutionFailedResponse])

    // Super-helpful debug in case the metadata watcher is still unhappy:
    if(metadataWatcherActor.underlyingActor.unsatisfiedMatchers.nonEmpty) {
      requiredMetadataMatchers foreach { matcher =>
        matcher.nearMissInformation.foreach { info =>
          System.out.println("A matcher had a near miss (it might still get a matching value later!): " + info.replace("\n", "..."))
        }
      }
    }

    // TODO: Yes, this might be slow... I'd advocate for refactoring away from the run-a-wdl style, but (shrug)
    // (but in fact, this never really takes 2 minutes. That's just for safety)
    Await.result(awaitable = metadataSuccessPromise.future, atMost = 2.minutes.dilated)

    system.stop(serviceRegistryActor)
  }

  it should "execute a workflow with scatters" in {
    val serviceRegistry = mockServiceRegistryActor
    val jobStore = system.actorOf(AlwaysHappyJobStoreActor.props)
    val subWorkflowStoreActor = system.actorOf(AlwaysHappySubWorkflowStoreActor.props)
    val callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props)
    val callCacheWriteActor = system.actorOf(EmptyCallCacheWriteActor.props)
    val dockerHashActor = system.actorOf(EmptyDockerHashActor.props)
    val ioActor = system.actorOf(SimpleIoActor.props)
    val jobTokenDispenserActor = system.actorOf(JobExecutionTokenDispenserActor.props)

    val MockBackendConfigEntry = BackendConfigurationEntry(
      name = MockBackendName,
      lifecycleActorFactoryClass = "cromwell.engine.backend.mock.DefaultBackendLifecycleActorFactory",
      stubbedConfig
    )
    CromwellBackends.initBackends(List(MockBackendConfigEntry))

    val workflowId = WorkflowId.randomId()
    val engineWorkflowDescriptor = createMaterializedEngineWorkflowDescriptor(workflowId, SampleWdl.SimpleScatterWdl.asWorkflowSources(runtime = runtimeSection))
    val workflowExecutionActor = system.actorOf(
      WorkflowExecutionActor.props(engineWorkflowDescriptor, ioActor, serviceRegistry, jobStore, subWorkflowStoreActor,
        callCacheReadActor, callCacheWriteActor, dockerHashActor, jobTokenDispenserActor, MockBackendSingletonCollection, AllBackendInitializationData.empty, restarting = false),
      "WorkflowExecutionActor")

    val scatterLog = "Starting calls: scatter0.inside_scatter:0:1, scatter0.inside_scatter:1:1, scatter0.inside_scatter:2:1, scatter0.inside_scatter:3:1, scatter0.inside_scatter:4:1"

    EventFilter.info(pattern = ".*Final Outputs", occurrences = 1).intercept {
      EventFilter.info(pattern = scatterLog, occurrences = 1).intercept {
        EventFilter.info(pattern = "Starting calls: scatter0.outside_scatter:NA:1", occurrences = 1).intercept {
          workflowExecutionActor ! ExecuteWorkflowCommand
        }
      }
    }
    system.stop(serviceRegistry)
  }
}
