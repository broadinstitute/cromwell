package cromwell.engine.backend

import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowId, WorkflowSourceFiles, WorkflowDescriptor}
import cromwell.engine.backend.local.LocalBackend
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpec, Matchers}

class BackendCallSpec extends CromwellTestkitSpec("BackendCallSpec") {

  val backend = new LocalBackend(new CromwellTestkitSpec.TestWorkflowManagerSystem().actorSystem)
  val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources()
  val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), sources)
  val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
  val backendCall = backend.bindCall(descriptor, CallKey(call, None), descriptor.actualInputs, AbortRegistrationFunction(_ => ()))

  "BackendCall hash function" should {
    "not change very often" in {
      val actual = backendCall.hash
      val expected = "20e7ec16a9d743a927791d530cba1807"
      assert(actual == expected,  s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
    }
  }
}
