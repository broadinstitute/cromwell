package cromwell.engine.backend

import cromwell.CallCachingWorkflowSpec._
import cromwell.CromwellTestkitSpec
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.workflow.BackendCallKey
import cromwell.util.{CromwellAggregatedException, SampleWdl}
import org.scalatest.concurrent.ScalaFutures
import wdl4s.values.WdlFile

class BackendCallSpec extends CromwellTestkitSpec with ScalaFutures with WorkflowDescriptorBuilder {
  override implicit val actorSystem = system

  val backend = new LocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, new CromwellTestkitSpec.TestWorkflowManagerSystem().actorSystem)

  "BackendCall hash function" should {
    "not hash by default" in {
      val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources()
      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources)
      val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
      val jobDescriptor = BackendCallJobDescriptor(descriptor, BackendCallKey(call, None, 1), descriptor.actualInputs)
      val actual = jobDescriptor.hash.futureValue.overallHash
      actual should be(empty)
    }

    "not change very often - if it changes, make sure it is for a good reason" in {
      val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources()
      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources, conf = callCachingConfig)
      val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
      val jobDescriptor = BackendCallJobDescriptor(descriptor, BackendCallKey(call, None, 1), descriptor.actualInputs)
      val actual = jobDescriptor.hash.futureValue.overallHash
      val expected = "1e5fc2d1fb3c8a26add14c3a5813b507"
      assert(actual == expected, s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
    }

    "not change with a Docker image hash specified - if it changes, make sure it is for a good reason" in {
      val nameAndDigest = "ubuntu@sha256:a2c950138e95bf603d919d0f74bec16a81d5cc1e3c3d574e8d5ed59795824f47"
      val sources = SampleWdl.CallCachingHashingWdl.asWorkflowSources( s"""runtime { docker: "$nameAndDigest" } """)
      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources, conf = callCachingConfig)
      val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "t").get
      val jobDescriptor = BackendCallJobDescriptor(descriptor, BackendCallKey(call, None, 1), descriptor.actualInputs)

      val actual = jobDescriptor.hash.futureValue.overallHash
      val expected = "e0fac3e3f6d3b611334fd2e0a504e99a"
      assert(actual == expected, s"Expected BackendCall hash to be $expected, but got $actual.  Did the hashing algorithm change?")
    }
  }

  "BackendCall adjustSharedInputPaths" should {
    "produce the correctly formatted exception message" in {
      val sources = SampleWdl.CannedFilePassing.asWorkflowSources()
      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = sources)
      val call = descriptor.namespace.workflow.calls.find(_.unqualifiedName == "wc").get
      val inputs = Map("in_file" -> WdlFile("/bad/path/to/file"))
      val jobDescriptor = BackendCallJobDescriptor(descriptor, BackendCallKey(call, None, 1), inputs)
      val exception = intercept[CromwellAggregatedException](backend.adjustSharedInputPaths(jobDescriptor))
      exception.getMessage should startWith("Failures during localization\nCould not localize /bad/path/to/file -> ")
      exception.getMessage should endWith("/call-wc/bad/path/to/file")
    }
  }
}
