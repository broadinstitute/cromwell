package cromwell.engine.backend

import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor, WorkflowId}
import cromwell.util.SampleWdl
import org.scalatest.concurrent.ScalaFutures

import scala.concurrent.ExecutionContext.Implicits.global


//TODO: This file should be removed imho since we do not plan to have backendCall
//class BackendCallSpec extends CromwellTestkitSpec with ScalaFutures {
//
//  val backend = new LocalBackend()
//  val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources()
//  val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), sources)
//  val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
//  val backendCall = backend.bindCall(descriptor, CallKey(call, None), descriptor.actualInputs, AbortRegistrationFunction(_ => ()))
//
//  "BackendCall hash function" should {
//    "not change very often - if it changes, make sure it is for a good reason" in {
//      val actual = backendCall.hash.futureValue.overallHash
//      val expected = "fe71298bb881a586178dca7c92fa945f"
//      assert(actual == expected, s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
//    }
//
//    "not change with a Docker image hash specified - if it changes, make sure it is for a good reason" in {
//      val nameAndDigest = "ubuntu@sha256:a2c950138e95bf603d919d0f74bec16a81d5cc1e3c3d574e8d5ed59795824f47"
//      val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources( s"""runtime { docker: "$nameAndDigest" } """)
//      val descriptor = WorkflowDescriptor(WorkflowId(UUID.randomUUID()), sources)
//      val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
//      val backendCall = backend.bindCall(descriptor, CallKey(call, None), descriptor.actualInputs, AbortRegistrationFunction(_ => ()))
//
//      val actual = backendCall.hash.futureValue.overallHash
//      val expected = "ca6ee457780b78290785f112c3a3acb4"
//      assert(actual == expected, s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
//    }
//  }
//}
