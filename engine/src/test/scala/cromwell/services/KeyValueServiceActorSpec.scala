package cromwell.services

import akka.pattern.ask
import com.typesafe.config.ConfigFactory
import cromwell.engine.{WorkflowSucceeded, WorkflowSourceFiles}
import cromwell.services.KeyValueServiceActor._

class KeyValueServiceActorSpec extends CromwellServicesSpec {

  val cromwellConfig = ConfigFactory.parseString(
    s"""services: {
        |  KeyValue: {
        |    class: "cromwell.services.KeyValueServiceActor"
        |    config {
        |      option1: "value1"
        |    }
        |  }
        |}
     """.stripMargin
  )

  val sources = WorkflowSourceFiles(
    wdlSource="""task a {command{}}
                |workflow w {
                |  call a
                |  call a as b
                |}
              """.stripMargin,
    inputsJson="{}",
    workflowOptionsJson="{}"
  )

  it should "insert a key/value" in {
    val descriptor = makeWorkflowDescriptor(sources)
    val kvActor = makeKeyValueActor(cromwellConfig, descriptor)
    val callKey = getBackendJobDescriptorKey(descriptor, "a")
    val put = KvPut(KvPair(ScopedKey(callKey, "k"), Option("v")))
    val get = KvGet(ScopedKey(callKey, "k"))

    (for {
      _ <- dataAccess.createWorkflow(descriptor, Nil, descriptor.namespace.workflow.calls, localBackend)
      put0 <- ask(kvActor, put).mapTo[KvResponse]
      _ = put0 shouldEqual KvPutSuccess(put)
      get0 <- ask(kvActor, get).mapTo[KvResponse]
      _ = get0 shouldEqual KvPair(ScopedKey(callKey, "k"), Option("v"))
      _ <- dataAccess.updateWorkflowState(descriptor.id, WorkflowSucceeded)
    } yield ()).futureValue
  }

  it should "return error if key doesn't exist" in {
    val descriptor = makeWorkflowDescriptor(sources)
    val kvActor = makeKeyValueActor(cromwellConfig, descriptor)
    val callKey = getBackendJobDescriptorKey(descriptor, "a")
    val get = KvGet(ScopedKey(callKey, "k"))

    (for {
      _ <- dataAccess.createWorkflow(descriptor, Nil, descriptor.namespace.workflow.calls, localBackend)
      get0 <- ask(kvActor, get).mapTo[KvKeyLookupFailed]
      _ = get0.action shouldEqual get
      _ <- dataAccess.updateWorkflowState(descriptor.id, WorkflowSucceeded)
    } yield ()).futureValue
  }

  it should "be able to overwrite values" in {
    val descriptor = makeWorkflowDescriptor(sources)
    val kvActor = makeKeyValueActor(cromwellConfig, descriptor)
    val callAKey = getBackendJobDescriptorKey(descriptor, "a")
    val putA = KvPut(KvPair(ScopedKey(callAKey, "k"), None))
    val putB = KvPut(KvPair(ScopedKey(callAKey, "k"), Option("v")))
    val get = KvGet(ScopedKey(callAKey, "k"))

    (for {
      _ <- dataAccess.createWorkflow(descriptor, Nil, descriptor.namespace.workflow.calls, localBackend)
      put0 <- ask(kvActor, putA).mapTo[KvPutSuccess]
      _ = put0 shouldEqual KvPutSuccess(putA)
      get0 <- ask(kvActor, get).mapTo[KvPair]
      _ = get0 shouldEqual KvPair(ScopedKey(callAKey, "k"), None)
      put1 <- ask(kvActor, putB).mapTo[KvPutSuccess]
      _ = put1 shouldEqual KvPutSuccess(putB)
      get1 <- ask(kvActor, get).mapTo[KvPair]
      _ = get1 shouldEqual KvPair(ScopedKey(callAKey, "k"), Option("v"))
      _ <- dataAccess.updateWorkflowState(descriptor.id, WorkflowSucceeded)
    } yield ()).futureValue
  }

  it should "be able to store NULL values" in {
    val descriptor = makeWorkflowDescriptor(sources)
    val kvActor = makeKeyValueActor(cromwellConfig, descriptor)
    val callAKey = getBackendJobDescriptorKey(descriptor, "a")
    val put = KvPut(KvPair(ScopedKey(callAKey, "k"), None))
    val get = KvGet(ScopedKey(callAKey, "k"))

    (for {
      _ <- dataAccess.createWorkflow(descriptor, Nil, descriptor.namespace.workflow.calls, localBackend)
      put0 <- ask(kvActor, put).mapTo[KvPutSuccess]
      _ = put0 shouldEqual KvPutSuccess(put)
      get0 <- ask(kvActor, get).mapTo[KvPair]
      _ = get0 shouldEqual KvPair(ScopedKey(callAKey, "k"), None)
      _ <- dataAccess.updateWorkflowState(descriptor.id, WorkflowSucceeded)
    } yield ()).futureValue
  }

  it should "partition keys by call" in {
    val descriptor = makeWorkflowDescriptor(sources)
    val kvActor = makeKeyValueActor(cromwellConfig, descriptor)
    val callAKey = getBackendJobDescriptorKey(descriptor, "a")
    val callBKey = getBackendJobDescriptorKey(descriptor, "b")
    val putA = KvPut(KvPair(ScopedKey(callAKey, "k"), Option("value_a")))
    val putB = KvPut(KvPair(ScopedKey(callBKey, "k"), Option("value_b")))
    val getA = KvGet(ScopedKey(callAKey, "k"))
    val getB = KvGet(ScopedKey(callBKey, "k"))

    (for {
      _ <- dataAccess.createWorkflow(descriptor, Nil, descriptor.namespace.workflow.calls, localBackend)
      put0 <- ask(kvActor, putA).mapTo[KvPutSuccess]
      _ = put0 shouldEqual KvPutSuccess(putA)
      get0 <- ask(kvActor, getB).mapTo[KvKeyLookupFailed]
      _ = get0.action shouldEqual getB
      put1 <- ask(kvActor, putB).mapTo[KvPutSuccess]
      _ = put1 shouldEqual KvPutSuccess(putB)
      get1 <- ask(kvActor, getA).mapTo[KvPair]
      _ = get1 shouldEqual KvPair(ScopedKey(callAKey, "k"), Option("value_a"))
      get2 <- ask(kvActor, getB).mapTo[KvPair]
      _ = get2 shouldEqual KvPair(ScopedKey(callBKey, "k"), Option("value_b"))
      _ <- dataAccess.updateWorkflowState(descriptor.id, WorkflowSucceeded)
    } yield ()).futureValue
  }

  it should "partition keys by workflow" in {
    val descriptor0 = makeWorkflowDescriptor(sources)
    val kv0 = makeKeyValueActor(cromwellConfig, descriptor0)
    val descriptor1 = makeWorkflowDescriptor(sources)
    val kv1 = makeKeyValueActor(cromwellConfig, descriptor1)

    val callA0 = getBackendJobDescriptorKey(descriptor0, "a")
    val callA1 = getBackendJobDescriptorKey(descriptor1, "a")
    val putA0 = KvPut(KvPair(ScopedKey(callA0, "k"), Option("v")))
    val putA1 = KvPut(KvPair(ScopedKey(callA1, "k"), Option("v")))
    val getA0 = KvGet(ScopedKey(callA0, "k"))
    val getA1 = KvGet(ScopedKey(callA1, "k"))

    (for {
      _ <- dataAccess.createWorkflow(descriptor0, Nil, descriptor0.namespace.workflow.calls, localBackend)
      _ <- dataAccess.createWorkflow(descriptor1, Nil, descriptor1.namespace.workflow.calls, localBackend)
      put0 <- ask(kv0, putA0).mapTo[KvPutSuccess]
      _ = put0 shouldEqual KvPutSuccess(putA0)
      get0 <- ask(kv1, getA1).mapTo[KvKeyLookupFailed]
      _ = get0.action shouldEqual getA1
      put1 <- ask(kv1, putA1).mapTo[KvPutSuccess]
      _ = put1 shouldEqual KvPutSuccess(putA1)
      get1 <- ask(kv0, getA0).mapTo[KvPair]
      _ = get1 shouldEqual KvPair(ScopedKey(callA0, "k"), Option("v"))
      get2 <- ask(kv1, getA1).mapTo[KvPair]
      _ = get2 shouldEqual KvPair(ScopedKey(callA1, "k"), Option("v"))
      _ <- dataAccess.updateWorkflowState(descriptor0.id, WorkflowSucceeded)
      _ <- dataAccess.updateWorkflowState(descriptor1.id, WorkflowSucceeded)
    } yield ()).futureValue
  }
}