package cromwell.engine.backend

import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{WorkflowDescriptor, WorkflowId}
import cromwell.util.SampleWdl
import org.scalatest.concurrent.ScalaFutures


class BackendCallSpec extends CromwellTestkitSpec with ScalaFutures {

  val backend = new LocalBackend(new CromwellTestkitSpec.TestWorkflowManagerSystem().actorSystem)
  val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources()
  val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), sources)
  val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
  val backendCall = backend.bindCall(descriptor, BackendCallKey(call, None), descriptor.actualInputs, abortRegistrationFunction = None)

  "BackendCall hash function" should {
    "not change very often - if it changes, make sure it is for a good reason" in {
      val actual = backendCall.hash.futureValue.overallHash
      val expected = "f91d97758ef2a97e9a4af15ad61bf1e8"
      assert(actual == expected, s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
    }

    "not change with a Docker image hash specified - if it changes, make sure it is for a good reason" in {
      val nameAndDigest = "ubuntu@sha256:a2c950138e95bf603d919d0f74bec16a81d5cc1e3c3d574e8d5ed59795824f47"
      val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources( s"""runtime { docker: "$nameAndDigest" } """)
      val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), sources)
      val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
      val backendCall = backend.bindCall(descriptor, BackendCallKey(call, None), descriptor.actualInputs, abortRegistrationFunction = None)

      val actual = backendCall.hash.futureValue.overallHash
      val expected = "711b81df90e257420ff571ac05bfdabf"
      assert(actual == expected, s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
    }
  }
}
