package cromwell.services

import akka.pattern.ask
import com.typesafe.config.ConfigFactory
import cromwell.engine.{WorkflowSourceFiles, WorkflowSucceeded}
import cromwell.services.KeyValueServiceActor._
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

class ServiceRegistryActorSpec extends FlatSpec with Matchers with ScalaFutures with CromwellServicesSpec {

  val cromwellConfig = ConfigFactory.parseString(
    s"""services: {
        |  KeyValue: {
        |    class: "cromwell.services.KeyValueServiceActor"
        |  }
        |}
     """.stripMargin
  )

  val sources = WorkflowSourceFiles(
    wdlSource="task a {command{}} workflow w {call a}",
    inputsJson="{}",
    workflowOptionsJson="{}"
  )

  it should "insert a key/value" in {
    val serviceRegistryActor = makeServiceRegistryActor(dataAccess, cromwellConfig)
    val descriptor = makeWorkflowDescriptor(sources)
    val callKey = getBackendJobDescriptorKey(descriptor, "a")
    val kvPut = KvPut(KvPair(ScopedKey(callKey, "k"), Option("v")))
    val kvGet = KvGet(ScopedKey(callKey, "k"))

    val future = for {
      _ <- dataAccess.createWorkflow(descriptor, Nil, descriptor.namespace.workflow.calls, localBackend)
      get0 <- ask(serviceRegistryActor, kvGet).mapTo[KvKeyLookupFailed]
      _ = get0.action shouldEqual kvGet
      put0 <- ask(serviceRegistryActor, kvPut).mapTo[KvResponse]
      _ = put0 shouldEqual KvPutSuccess(kvPut)
      get1 <- ask(serviceRegistryActor, kvGet).mapTo[KvResponse]
      _ = get1 shouldEqual KvPair(ScopedKey(callKey, "k"), Option("v"))
      // Put the workflow in a terminal state so tests running after it don't restart it and become confused.
      _ <- dataAccess.updateWorkflowState(descriptor.id, WorkflowSucceeded)
    } yield ()

    future.futureValue
  }
}