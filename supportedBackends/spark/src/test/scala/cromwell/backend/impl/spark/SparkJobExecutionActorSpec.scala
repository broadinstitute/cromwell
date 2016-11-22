package cromwell.backend.impl.spark

import java.io.Writer
import java.nio.file.Path

import akka.testkit.{ImplicitSender, TestActorRef}
import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{JobFailedNonRetryableResponse, JobSucceededResponse}
import cromwell.backend.impl.spark.SparkClusterProcess._
import cromwell.backend.io._
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendSpec}
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.core.path.{PathWriter, TailedWriter, UntailedWriter}
import org.mockito.Mockito
import org.mockito.Mockito._
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.mockito.MockitoSugar
import org.scalatest.{BeforeAndAfter, Matchers, WordSpecLike}
import wdl4s._
import wdl4s.values.WdlValue

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.sys.process.{Process, ProcessLogger}

class SparkJobExecutionActorSpec extends TestKitSuite("SparkJobExecutionActor")
  with WordSpecLike
  with Matchers
  with MockitoSugar
  with BeforeAndAfter
  with ImplicitSender {

  import BackendSpec._

  private val sparkProcess: SparkProcess = mock[SparkProcess]
  private val sparkCommands: SparkCommands = new SparkCommands
  private val sparkClusterProcess: SparkClusterProcess = mock[SparkClusterProcess]

  private val helloWorldWdl =
    """
      |task hello {
      |
      |  command {
      |    sparkApp
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |  RUNTIME
      |}
      |
      |workflow wf_hello {
      |  call hello
      |}
    """.stripMargin

  private val helloWorldClusterWdl =
    """
      |task helloClusterMode {
      |
      | command {
      |   sparkApp
      |  }
      | output {
      |   String salutation = read_string(stdout())
      | }
      | RUNTIME
      |}
      |
      |workflow wf_helloClusterMode {
      |   call helloClusterMode
      |}
    """.stripMargin

  private val sampleSubmissionResponse =
    """
      |Running Spark using the REST application submission protocol.
      |16/08/06 18:35:26 INFO rest.RestSubmissionClient: Submitting a request to launch an application in spark://host-10-0-1-53:6066.
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: Submission successfully created as driver-20160806183527-0006. Polling submission state...
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: Submitting a request for the status of submission driver-20160806183527-0006 in spark://host-10-0-1-53:6066.
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: State of driver driver-20160806183527-0006 is now RUNNING.
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: Driver is running on worker worker-20160801162431-10.0.1.55-43834 at 10.0.1.55:43834.
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: Server responded with CreateSubmissionResponse:
      |{
      |  "action" : "CreateSubmissionResponse",
      |  "message" : "Driver successfully submitted as driver-20160806183527-0006",
      |  "serverSparkVersion" : "1.6.1",
      |  "submissionId" : "driver-20160806183527-0006",
      |  "success" : true
      |}
    """.stripMargin

  private val failedOnStderr =
    """
      |runtime {
      | appMainClass: "test"
      | failOnStderr: true
      |}
    """.stripMargin

  private val passOnStderr =
    """
      |runtime {
      | appMainClass: "test"
      | failOnStderr: false
      |}
    """.stripMargin

  private val backendClusterConfig = ConfigFactory.parseString(
    s"""{
        |  root = "local-cromwell-executions"
        |  filesystems {
        |    local {
        |      localization = [
        |        "hard-link", "soft-link", "copy"
        |      ]
        |    }
        |  }
        |  master: "spark"
        |  deployMode: "cluster"
        |}
        """.stripMargin)

  private val backendClientConfig = ConfigFactory.parseString(
    s"""{
        |  root = "local-cromwell-executions"
        |  filesystems {
        |    local {
        |      localization = [
        |        "hard-link", "soft-link", "copy"
        |      ]
        |    }
        |  }
        |  master: "local"
        |  deployMode: "client"
        |}
        """.stripMargin)

  private val timeout = Timeout(20.seconds)

  after {
    Mockito.reset(sparkProcess)
    Mockito.reset(sparkClusterProcess)
  }

  override def afterAll(): Unit = {
    system.terminate()
    ()
  }

  "executeTask method in cluster deploy mode " should {
    "return succeed response when the spark cluster process monitor method returns finished status" in {
      val jobDescriptor = prepareJob(helloWorldClusterWdl, passOnStderr, isCluster = true)
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderr = ""

      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val clusterExtProcess = sparkClusterProcess
        override lazy val cmds = sparkCommands
      }).underlyingActor

      when(sparkClusterProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(sparkClusterProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkClusterProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkClusterProcess.processStderr).thenReturn(stderr)
      when(sparkClusterProcess.startMonitoringSparkClusterJob(any[Path], any[String])).thenReturn(Future.successful(Finished))

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobSucceededResponse]
        verify(sparkClusterProcess, times(1)).externalProcess(any[Seq[String]], any[ProcessLogger])
        verify(sparkClusterProcess, times(1)).tailedWriter(any[Int], any[Path])
        verify(sparkClusterProcess, times(1)).untailedWriter(any[Path])
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(any[Path], any[String])
      }
      cleanUpJob(jobPaths)
    }

    "return failed non-retryable response when the spark monitor function returns failed status" in {
      val jobDescriptor = prepareJob(helloWorldClusterWdl, passOnStderr, isCluster = true)
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderr = ""

      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val clusterExtProcess = sparkClusterProcess
        override lazy val cmds = sparkCommands
      }).underlyingActor

      when(sparkClusterProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(sparkClusterProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkClusterProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkClusterProcess.processStderr).thenReturn(stderr)
      when(sparkClusterProcess.startMonitoringSparkClusterJob(any[Path], any[String])).thenReturn(Future.successful(Failed(new Throwable("failed to monitor"))))

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains("failed to monitor"))
        verify(sparkClusterProcess, times(1)).externalProcess(any[Seq[String]], any[ProcessLogger])
        verify(sparkClusterProcess, times(1)).tailedWriter(any[Int], any[Path])
        verify(sparkClusterProcess, times(1)).untailedWriter(any[Path])
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(any[Path], any[String])
      }
      cleanUpJob(jobPaths)
    }

    "return failed non-retryable response when the spark monitor function throws an exception" in {
      val jobDescriptor = prepareJob(helloWorldClusterWdl, passOnStderr, isCluster = true)
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderr = ""

      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val clusterExtProcess = sparkClusterProcess
        override lazy val cmds = sparkCommands
      }).underlyingActor

      when(sparkClusterProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(sparkClusterProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkClusterProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkClusterProcess.processStderr).thenReturn(stderr)
      when(sparkClusterProcess.startMonitoringSparkClusterJob(any[Path], any[String])).thenReturn(Future.failed(new IllegalStateException("failed to start monitoring process")))

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains("failed to start monitoring process"))
        verify(sparkClusterProcess, times(1)).externalProcess(any[Seq[String]], any[ProcessLogger])
        verify(sparkClusterProcess, times(1)).tailedWriter(any[Int], any[Path])
        verify(sparkClusterProcess, times(1)).untailedWriter(any[Path])
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(any[Path], any[String])
      }
      cleanUpJob(jobPaths)
    }

    "return failed non-retryable response when there's a zero exit code but stderr is  not empty" in {
      val jobDescriptor = prepareJob(helloWorldClusterWdl, failedOnStderr, isCluster = true)
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      File(jobPaths.stderr) < "failed"
      File(jobPaths.callExecutionRoot.resolve("cluster.json")) < sampleSubmissionResponse

      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val clusterExtProcess = sparkClusterProcess
        override lazy val cmds = sparkCommands
      }).underlyingActor

      when(sparkClusterProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(sparkClusterProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkClusterProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkClusterProcess.processStderr).thenReturn(sampleSubmissionResponse)

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains(s"Execution process failed although return code is zero but stderr is not empty"))
        verify(sparkClusterProcess, times(1)).externalProcess(any[Seq[String]], any[ProcessLogger])
        verify(sparkClusterProcess, times(1)).tailedWriter(any[Int], any[Path])
        verify(sparkClusterProcess, times(1)).untailedWriter(any[Path])
      }
      cleanUpJob(jobPaths)
    }

    "return failed non-retryable response when non zero exit code is received for submitting job" in {
      val jobDescriptor = prepareJob(helloWorldClusterWdl, passOnStderr, isCluster = true)
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""
      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val clusterExtProcess = sparkClusterProcess
        override lazy val cmds = sparkCommands
      }).underlyingActor

      when(sparkClusterProcess.commandList(any[String])).thenReturn(Seq.empty[String])
      when(sparkClusterProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(-1)
      when(sparkClusterProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkClusterProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkClusterProcess.processStderr).thenReturn(stderrResult)

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains(s"Execution process failed. Spark returned non zero status code:"))
      }
      cleanUpJob(jobPaths)
    }

    "return failed non-retryable response when submit job process exit with failure response" in {
      val jobDescriptor = prepareJob(helloWorldClusterWdl, passOnStderr, isCluster = true)
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""
      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val clusterExtProcess = sparkClusterProcess
        override lazy val cmds = sparkCommands
      }).underlyingActor

      when(sparkClusterProcess.commandList(any[String])).thenReturn(Seq.empty[String])
      when(sparkClusterProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenThrow(new IllegalStateException("submit job process exitValue method failed"))
      when(sparkClusterProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkClusterProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkClusterProcess.processStderr).thenReturn(stderrResult)

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains(s"submit job process exitValue method failed"))
      }
      cleanUpJob(jobPaths)
    }
  }

  "executeTask method in Spark client and Yarn client, cluster mode" should {
    "return succeeded task status with stdout" in {
      val jobDescriptor = prepareJob()
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""

      when(sparkProcess.commandList(any[String])).thenReturn(Seq.empty[String])
      when(sparkProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(sparkProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkProcess.processStderr).thenReturn(stderrResult)

      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val extProcess = sparkProcess
        override lazy val cmds = sparkCommands
      }).underlyingActor

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobSucceededResponse]
        verify(sparkProcess, times(1)).externalProcess(any[Seq[String]], any[ProcessLogger])
        verify(sparkProcess, times(1)).tailedWriter(any[Int], any[Path])
        verify(sparkProcess, times(1)).untailedWriter(any[Path])
      }

      cleanUpJob(jobPaths)
    }

    "return failed task status on non-zero process exit" in {
      val jobDescriptor = prepareJob()
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)

      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val cmds = sparkCommands
        override lazy val extProcess = sparkProcess
      }).underlyingActor
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""

      when(sparkProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(-1)
      when(sparkProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkProcess.processStderr).thenReturn(stderrResult)

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains(s"Execution process failed. Spark returned non zero status code:"))
      }

      cleanUpJob(jobPaths)
    }

    "return failed task status when stderr is non-empty but process exit with zero return code" in {
      val jobDescriptor = prepareJob(runtimeString = failedOnStderr)
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)

      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val cmds = sparkCommands
        override lazy val extProcess = sparkProcess
      }).underlyingActor
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      File(jobPaths.stderr) < "failed"

      when(sparkProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(sparkProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains(s"Execution process failed although return code is zero but stderr is not empty"))
      }

      cleanUpJob(jobPaths)
    }

    "return succeeded task status when stderr is non-empty but process exit with zero return code" in {
      val jobDescriptor = prepareJob()
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)

      val backend = TestActorRef(new SparkJobExecutionActor(job, backendConfigDesc) {
        override lazy val cmds = sparkCommands
        override lazy val extProcess = sparkProcess
      }).underlyingActor
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter

      when(sparkProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(sparkProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobSucceededResponse]
        verify(sparkProcess, times(1)).externalProcess(any[Seq[String]], any[ProcessLogger])
        verify(sparkProcess, times(1)).tailedWriter(any[Int], any[Path])
        verify(sparkProcess, times(1)).untailedWriter(any[Path])
      }

      cleanUpJob(jobPaths)
    }

  }

  private def cleanUpJob(jobPaths: JobPathsWithDocker): Unit = {
    File(jobPaths.workflowRoot).delete(true)
    ()
  }

  private def prepareJob(wdlSource: WdlSource = helloWorldWdl, runtimeString: String = passOnStderr, inputFiles: Option[Map[String, WdlValue]] = None, isCluster: Boolean = false): TestJobDescriptor = {
    val backendWorkflowDescriptor = buildWorkflowDescriptor(wdl = wdlSource, inputs = inputFiles.getOrElse(Map.empty), runtime = runtimeString)
    val backendConfigurationDescriptor = if (isCluster) BackendConfigurationDescriptor(backendClusterConfig, ConfigFactory.load) else BackendConfigurationDescriptor(backendClientConfig, ConfigFactory.load)
    val jobDesc = jobDescriptorFromSingleCallWorkflow(backendWorkflowDescriptor, inputFiles.getOrElse(Map.empty), WorkflowOptions.empty, Set.empty)
    val jobPaths = if (isCluster) new JobPathsWithDocker(jobDesc.key, backendWorkflowDescriptor, backendClusterConfig) else new JobPathsWithDocker(jobDesc.key, backendWorkflowDescriptor, backendClientConfig)
    val executionDir = jobPaths.callExecutionRoot
    val stdout = File(executionDir.toString, "stdout")
    stdout.createIfNotExists(asDirectory = false, createParents = true)
    val stderr = File(executionDir.toString, "stderr")
    stderr.createIfNotExists(asDirectory = false, createParents = true)
    TestJobDescriptor(jobDesc, jobPaths, backendConfigurationDescriptor)
  }

  private case class TestJobDescriptor(jobDescriptor: BackendJobDescriptor, jobPaths: JobPathsWithDocker, backendConfigurationDescriptor: BackendConfigurationDescriptor)

  trait MockWriter extends Writer {
    var closed = false

    override def close() = closed = true

    override def flush() = {}

    override def write(a: Array[Char], b: Int, c: Int) = {}
  }

  trait MockPathWriter extends PathWriter {
    override lazy val writer: Writer = new MockWriter {}
    override val path: Path = mock[Path]
  }

}
