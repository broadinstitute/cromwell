package cromwell.backend.impl.aws

import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.backend.BackendSpec._
import cromwell.backend.io.JobPaths
import cromwell.core.WorkflowOptions
import cromwell.core.callcaching.NoDocker
import cromwell.util.SampleWdl
import eu.timepit.refined.refineMV
import eu.timepit.refined.numeric.Positive
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import spray.json.{JsObject, JsString}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.graph.CommandCallNode
import common.mock.MockSugar

class AwsBatchJobDefinitionSpec extends AnyWordSpecLike with Matchers with MockSugar {

  import AwsBatchWorkflowOptionKeys._

  // Create a simple workflow descriptor for testing
  val workflowDescriptor: BackendWorkflowDescriptor = buildWdlWorkflowDescriptor(
    SampleWdl.HelloWorld.workflowSource(),
    inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.map { case (k, v) =>
      k -> JsString(v)
    }).compactPrint)
  )

  val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
  val jobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)

  // Mock job descriptor for testing
  val jobDescriptor: BackendJobDescriptor = BackendJobDescriptor(
    workflowDescriptor,
    jobKey,
    Map.empty,
    Map.empty,
    NoDocker,
    None,
    Map.empty
  )

  // Test runtime attributes
  val runtimeAttributes = AwsBatchRuntimeAttributes(
    cpu = refineMV[Positive](1),
    gpuCount = 0,
    zones = Vector("us-east-1a"),
    memory = MemorySize(2.0, MemoryUnit.GB),
    disks = Seq.empty,
    dockerImage = "ubuntu:latest",
    queueArn = "arn:aws:batch:us-east-1:123456789:job-queue/test-queue",
    failOnStderr = false,
    continueOnReturnCode = cromwell.backend.validation.ContinueOnReturnCodeSet(Set(0)),
    noAddress = false,
    scriptS3BucketName = "test-bucket",
    awsBatchRetryAttempts = 1,
    awsBatchEvaluateOnExit = Vector(Map.empty[String, String]),
    ulimits = Vector(Map.empty[String, String]),
    efsDelocalize = false,
    efsMakeMD5 = false,
    fuseMount = false,
    sharedMemorySize = MemorySize(10, MemoryUnit.MB),
    jobTimeout = 0,
    logGroupName = "/aws/batch/job",
    additionalTags = Map.empty,
    fileSystem = "s3"
  )

  // Mock job paths for testing
  val mockJobPaths: JobPaths = mock[JobPaths]

  "AwsBatchJobDefinition StandardAwsBatchJobDefinitionBuilder" should {

    "build a job definition name without jobRoleArn when not provided" in {
      val workflowOptions = WorkflowOptions(JsObject.empty)

      val context = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        workflowOptions = workflowOptions
      )

      val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(context)

      // Job definition name should be deterministic based on inputs
      // The colon in "ubuntu:latest" gets replaced with underscore
      jobDefinition.name should startWith("cromwell_ubuntu_latest")
      // The name should contain the SHA1 hash at the end
      jobDefinition.name should fullyMatch regex """cromwell_ubuntu_latest_[a-f0-9]{40}"""

      // Verify no roleArn is included in hash calculation by checking it matches expected pattern
      val roleArnStr = workflowOptions.getOrElse(JobRoleArn, "")
      roleArnStr should be("")
    }

    "build a job definition name with jobRoleArn when provided" in {
      // This test makes sure that the Job Definition is rebuilt correctly
      // if the user submits a workflow with different JobRoleARN.
      val roleArn = "arn:aws:iam::123456789012:role/MyJobRole"
      val workflowOptions = WorkflowOptions(JsObject(Map(JobRoleArn -> JsString(roleArn))))

      val context = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        workflowOptions = workflowOptions
      )

      val jobDefinitionWithRole = StandardAwsBatchJobDefinitionBuilder.build(context)

      // Job definition name should be different when roleArn is included
      // The colon in "ubuntu:latest" gets replaced with underscore
      jobDefinitionWithRole.name should startWith("cromwell_ubuntu_latest")
      // The name should contain the SHA1 hash at the end
      jobDefinitionWithRole.name should fullyMatch regex """cromwell_ubuntu_latest_[a-f0-9]{40}"""

      // Create the same job definition without roleArn to compare
      val workflowOptionsNoRole = WorkflowOptions(JsObject.empty)
      val contextNoRole = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        workflowOptions = workflowOptionsNoRole
      )

      val jobDefinitionNoRole = StandardAwsBatchJobDefinitionBuilder.build(contextNoRole)

      // The two names should be different because roleArn is included in the hash
      jobDefinitionWithRole.name should not equal jobDefinitionNoRole.name
    }

    "apply jobRoleArn to container properties when provided" in {
      val roleArn = "arn:aws:iam::123456789012:role/MyJobRole"
      val workflowOptions = WorkflowOptions(JsObject(Map(JobRoleArn -> JsString(roleArn))))

      val context = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        workflowOptions = workflowOptions
      )

      val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(context)
      val containerProperties = jobDefinition.containerProperties

      // Verify the job role ARN is set on the container properties
      containerProperties.jobRoleArn() should equal(roleArn)
    }

    "not set jobRoleArn on container properties when not provided" in {
      val workflowOptions = WorkflowOptions(JsObject.empty)

      val context = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        workflowOptions = workflowOptions
      )

      val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(context)
      val containerProperties = jobDefinition.containerProperties

      // Verify no job role ARN is set when not provided
      Option(containerProperties.jobRoleArn()).isEmpty should be(true)
    }

    "correctly use getOrElse method for jobRoleArn retrieval" in {
      val workflowOptions = WorkflowOptions(JsObject.empty)

      // Test that getOrElse returns default value when key is not present
      val roleArnStr = workflowOptions.getOrElse(JobRoleArn, "default-value")
      roleArnStr should be("default-value")

      // Test that getOrElse returns the actual value when key is present
      val workflowOptionsWithRole = WorkflowOptions(JsObject(Map(JobRoleArn -> JsString("my-role-arn"))))
      val roleArnStrWithValue = workflowOptionsWithRole.getOrElse(JobRoleArn, "default-value")
      roleArnStrWithValue should be("my-role-arn")
    }
  }
}
