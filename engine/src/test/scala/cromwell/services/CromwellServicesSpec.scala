package cromwell.services

import akka.actor.ActorRef
import akka.util.Timeout
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.{WorkflowId, WorkflowOptions}
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.local.OldStyleLocalBackend
import cromwell.engine.backend.{OldStyleWorkflowDescriptor, WorkflowDescriptorBuilder}
import cromwell.engine.db.DataAccess
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import spray.json.JsObject

import scala.concurrent.duration._

trait CromwellServicesSpec extends FlatSpec with Matchers with BeforeAndAfterAll with WorkflowDescriptorBuilder with ScalaFutures {
  val workflowManagerSystem = new TestWorkflowManagerSystem
  val actorSystem = workflowManagerSystem.actorSystem
  val dataAccess = DataAccess.globalDataAccess

  val localBackend = OldStyleLocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, actorSystem)
  implicit val timeout = Timeout(5.seconds)
  implicit val ec = actorSystem.dispatcher
  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  override protected def afterAll() = {
    workflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  protected def getBackendJobDescriptorKey(descriptor: OldStyleWorkflowDescriptor, callName: String): BackendJobDescriptor = {
    val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == callName).get
    val backendWorkflowDescriptor = BackendWorkflowDescriptor(
      descriptor.id,
      descriptor.namespace,
      Map.empty,
      WorkflowOptions(JsObject())
    )
    val key = BackendJobDescriptorKey(call, None, 1)
    BackendJobDescriptor(backendWorkflowDescriptor, key, Map.empty)
  }

  protected def makeWorkflowDescriptor(sources: WorkflowSourceFiles): OldStyleWorkflowDescriptor = {
    val workflowId = WorkflowId.randomId()
    materializeWorkflowDescriptorFromSources(id = workflowId, workflowSources = sources)
  }

  protected def makeKeyValueActor(cromwellConfig: Config, workflowDescriptor: OldStyleWorkflowDescriptor): ActorRef = {
    actorSystem.actorOf(
      KeyValueServiceActor.props(ConfigFactory.parseString(""), cromwellConfig)
    )
  }

  protected def makeServiceRegistryActor(dataAccess: DataAccess, cromwellConfig: Config): ActorRef = {
    actorSystem.actorOf(
      ServiceRegistryActor.props(cromwellConfig)
    )
  }
}
