package cromwell.engine.workflow

import cromwell.engine.WorkflowDescriptor
import cromwell.engine.io.gcs.GcsFileSystem
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wdl4s.Call

import scala.util.{Failure, Success}

class BackendCallKeySpec extends FlatSpec with Matchers with Mockito {

  it should "build correct call root path" in {
    import cromwell.engine.PathString._

    import scala.concurrent.ExecutionContext.Implicits.global

    val call = mock[Call]
    call.unqualifiedName returns "callName"
    val wd = mock[WorkflowDescriptor]
    val root = "/root"
    val bucket = "gs://bucket"
    val gcsFileSystem = Success(GcsFileSystem(Failure(new Throwable("no interface")), bucket))
    wd.workflowRootPathWithBaseRoot(anyString) answers {root => s"$root/wfName/uuid".toPath(gcsFileSystem) }

    new BackendCallKey(call, None, 1).callRootPathWithBaseRoot(wd, root).toString shouldBe "/root/wfName/uuid/call-callName"
    new BackendCallKey(call, Some(0), 1).callRootPathWithBaseRoot(wd, root).toString shouldBe "/root/wfName/uuid/call-callName/shard-0"
    new BackendCallKey(call, Some(0), 2).callRootPathWithBaseRoot(wd, root).toString shouldBe "/root/wfName/uuid/call-callName/shard-0/attempt-2"
    new BackendCallKey(call, None, 2).callRootPathWithBaseRoot(wd, root).toString shouldBe "/root/wfName/uuid/call-callName/attempt-2"

    new BackendCallKey(call, None, 1).callRootPathWithBaseRoot(wd, bucket).toString shouldBe "gs://bucket/wfName/uuid/call-callName"
    new BackendCallKey(call, Some(0), 1).callRootPathWithBaseRoot(wd, bucket).toString shouldBe "gs://bucket/wfName/uuid/call-callName/shard-0"
    new BackendCallKey(call, Some(0), 2).callRootPathWithBaseRoot(wd, bucket).toString shouldBe "gs://bucket/wfName/uuid/call-callName/shard-0/attempt-2"
    new BackendCallKey(call, None, 2).callRootPathWithBaseRoot(wd, bucket).toString shouldBe "gs://bucket/wfName/uuid/call-callName/attempt-2"
  }

  it should "increment attempt attribute when cloning itself for retry" in {
    val call = mock[Call]
    new BackendCallKey(call, None, 1).retryClone.attempt shouldBe 2
  }
}
