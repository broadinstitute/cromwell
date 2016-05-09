package cromwell.engine.workflow

import cromwell.engine.backend.{OldStyleBackendCallJobDescriptor, OldStyleBackend, OldStyleWorkflowDescriptor}
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wdl4s.Call

class BackendCallKeySpec extends FlatSpec with Matchers with Mockito {

  it should "build correct call root path" in {
    import cromwell.engine.backend.io._

    val root = "/root"
    val bucket = "gs://bucket"
    val call = mock[Call]
    call.unqualifiedName returns "callName"
    val wd = mock[OldStyleWorkflowDescriptor]
    wd.fileSystems returns defaultFileSystems
    wd.workflowRootPathWithBaseRoot(anyString) answers { root => s"$root/wfName/uuid".toPath(defaultFileSystems) }

    OldStyleBackend.callRootPathWithBaseRoot(new OldStyleBackendCallJobDescriptor(
      wd, new BackendCallKey(call, None, 1)), root).toString shouldBe "/root/wfName/uuid/call-callName"
    OldStyleBackend.callRootPathWithBaseRoot(new OldStyleBackendCallJobDescriptor(
      wd, new BackendCallKey(call, Some(0), 1)), root).toString shouldBe "/root/wfName/uuid/call-callName/shard-0"
    OldStyleBackend.callRootPathWithBaseRoot(new OldStyleBackendCallJobDescriptor(
      wd, new BackendCallKey(call, Some(0), 2)), root).toString shouldBe "/root/wfName/uuid/call-callName/shard-0/attempt-2"
    OldStyleBackend.callRootPathWithBaseRoot(new OldStyleBackendCallJobDescriptor(
      wd, new BackendCallKey(call, None, 2)), root).toString shouldBe "/root/wfName/uuid/call-callName/attempt-2"
    OldStyleBackend.callRootPathWithBaseRoot(new OldStyleBackendCallJobDescriptor(
      wd, new BackendCallKey(call, None, 1)), bucket).toString shouldBe "gs://bucket/wfName/uuid/call-callName"
    OldStyleBackend.callRootPathWithBaseRoot(new OldStyleBackendCallJobDescriptor(
      wd, new BackendCallKey(call, Some(0), 1)), bucket).toString shouldBe "gs://bucket/wfName/uuid/call-callName/shard-0"
    OldStyleBackend.callRootPathWithBaseRoot(new OldStyleBackendCallJobDescriptor(
      wd, new BackendCallKey(call, Some(0), 2)), bucket).toString shouldBe "gs://bucket/wfName/uuid/call-callName/shard-0/attempt-2"
    OldStyleBackend.callRootPathWithBaseRoot(new OldStyleBackendCallJobDescriptor(
      wd, new BackendCallKey(call, None, 2)), bucket).toString shouldBe "gs://bucket/wfName/uuid/call-callName/attempt-2"
  }

  it should "increment attempt attribute when cloning itself for retry" in {
    val call = mock[Call]
    new BackendCallKey(call, None, 1).retryClone.attempt shouldBe 2
  }
}
