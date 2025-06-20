package cromwell.backend.google.batch.models

import com.google.cloud.NoCredentials
import common.collections.EnhancedCollections._
import cromwell.backend.BackendSpec
import cromwell.backend.google.batch.actors.GcpBatchInitializationActor
import cromwell.backend.google.batch.models.GcpBatchTestConfig.{
  batchAttributes,
  gcpBatchConfiguration,
  googleConfiguration,
  pathBuilders,
  BatchBackendConfigurationDescriptor
}
import cromwell.backend.io.JobPathsSpecHelper._
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import spray.json.{JsObject, JsString}

import scala.concurrent.ExecutionContext.Implicits.global

class GcpBatchJobPathsSpec extends TestKitSuite with AnyFlatSpecLike with Matchers {

  import BackendSpec._

  behavior of "GcpBatchCallPaths"

  it should "map the correct filenames" in {
    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)

    val workflowPaths = GcpBatchWorkflowPaths(
      workflowDescriptor,
      NoCredentials.getInstance(),
      NoCredentials.getInstance(),
      gcpBatchConfiguration,
      pathBuilders(),
      GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper
    )

    val callPaths = GcpBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.returnCodeFilename should be("rc")
    callPaths.stderr.getFileName.pathAsString should be("gs://my-cromwell-workflows-bucket/stderr")
    callPaths.stdout.getFileName.pathAsString should be("gs://my-cromwell-workflows-bucket/stdout")
    callPaths.batchLogFilename should be("hello.log")
  }

  it should "map the correct paths" in {
    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)

    val workflowPaths = GcpBatchWorkflowPaths(
      workflowDescriptor,
      NoCredentials.getInstance(),
      NoCredentials.getInstance(),
      gcpBatchConfiguration,
      pathBuilders(),
      GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper
    )

    val callPaths = GcpBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.returnCode.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/rc")
    callPaths.stdout.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stdout")
    callPaths.stderr.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stderr")
    callPaths.batchLogPath.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello.log")
  }

  it should "map the correct call context" in {

    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)

    val workflowPaths = GcpBatchWorkflowPaths(
      workflowDescriptor,
      NoCredentials.getInstance(),
      NoCredentials.getInstance(),
      gcpBatchConfiguration,
      pathBuilders(),
      GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper
    )

    val callPaths = GcpBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.callContext.root.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello")
    callPaths.callContext.stdout should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stdout")
    callPaths.callContext.stderr should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stderr")
  }

  it should "include the backend log path in metadata and detritus when the file is configured to exist" in {
    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)

    val configWithPathLogsPolicy =
      new GcpBatchConfiguration(BatchBackendConfigurationDescriptor,
                                googleConfiguration,
                                batchAttributes.copy(logsPolicy = GcpBatchLogsPolicy.Path)
      )

    val workflowPaths = GcpBatchWorkflowPaths(
      workflowDescriptor,
      NoCredentials.getInstance(),
      NoCredentials.getInstance(),
      configWithPathLogsPolicy,
      pathBuilders(),
      GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper
    )

    val callPaths = GcpBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.customMetadataPaths should contain(
      "backendLogs:log" -> callPaths.batchLogPath
    )
    callPaths.customDetritusPaths should contain(
      "log" -> callPaths.batchLogPath
    )
  }

  it should "not include the backend log path in metadata and detritus when the file is not configured to exist" in {
    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)

    val workflowPaths = GcpBatchWorkflowPaths(
      workflowDescriptor,
      NoCredentials.getInstance(),
      NoCredentials.getInstance(),
      gcpBatchConfiguration,
      pathBuilders(),
      GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper
    )

    val callPaths = GcpBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.customMetadataPaths shouldNot contain(
      "backendLogs:log" -> callPaths.batchLogPath
    )
    callPaths.customDetritusPaths shouldNot contain(
      "log" -> callPaths.batchLogPath
    )
  }
}
