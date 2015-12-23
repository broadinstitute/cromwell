package cromwell.engine.backend

import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor, WorkflowId}
import cromwell.util.SampleWdl
import org.scalatest.concurrent.ScalaFutures

import scala.concurrent.ExecutionContext.Implicits.global


class BackendCallSpec extends CromwellTestkitSpec with ScalaFutures {

  val backend = new LocalBackend(new CromwellTestkitSpec.TestWorkflowManagerSystem().actorSystem)
  val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources()
  val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), sources)
  val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
  val backendCall = backend.bindCall(descriptor, CallKey(call, None), descriptor.actualInputs, AbortRegistrationFunction(_ => ()))

  "BackendCall hash function" should {
    "not change very often" in {
      val actual = backendCall.hash.futureValue.overallHash
      val expected = "20e7ec16a9d743a927791d530cba1807"
      assert(actual == expected, s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
    }

    "not change with a Docker image hash specified" in {
      val nameAndDigest = "ubuntu@sha256:a2c950138e95bf603d919d0f74bec16a81d5cc1e3c3d574e8d5ed59795824f47"
      val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources( s"""runtime { docker: "$nameAndDigest" } """)
      val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), sources)
      val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
      val backendCall = backend.bindCall(descriptor, CallKey(call, None), descriptor.actualInputs, AbortRegistrationFunction(_ => ()))

      val actual = backendCall.hash.futureValue.overallHash
      val expected = "8ed446e4662cc7c0455231dc98713e87"
      assert(actual == expected, s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
    }
  }
}
