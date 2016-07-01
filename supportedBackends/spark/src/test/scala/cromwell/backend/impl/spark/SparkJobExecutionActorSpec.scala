package cromwell.backend.impl.spark

import java.io.{File, FileWriter, Writer}
import java.nio.file.{Files, Path, Paths}

import akka.actor.ActorSystem
import akka.testkit.{ImplicitSender, TestActorRef, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.backend.BackendJobExecutionActor.{FailedNonRetryableResponse, SucceededResponse}
import better.files._
import cromwell.backend.io.{BackendTestkitSpec, JobPaths}
import cromwell.core.{PathWriter, TailedWriter, UntailedWriter}
import org.mockito.Matchers._
import org.mockito.Mockito
import org.mockito.Mockito._
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.mock.MockitoSugar
import org.scalatest.{BeforeAndAfter, BeforeAndAfterAll, Matchers, WordSpecLike}
import wdl4s.values.{WdlFile, WdlValue}

import scala.concurrent.duration._
import scala.io.Source
import scala.sys.process.{Process, ProcessLogger}

class SparkJobExecutionActorSpec extends TestKit(ActorSystem("SparkJobExecutionActor"))
  with BackendTestkitSpec
  with WordSpecLike
  with Matchers
  with MockitoSugar
  with BeforeAndAfter
  with BeforeAndAfterAll
  with ImplicitSender {

  private val sparkProcess: SparkProcess = mock[SparkProcess]
  private val sparkCommands: SparkCommands = new SparkCommands

  private val helloWorldWdl =
    """
      |task hello {
      |
      |  command {
      |    echo "Hello World!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |  RUNTIME
      |}
      |
      |workflow hello {
      |  call hello
      |}
    """.stripMargin

  private val helloWorldWdlWithFileInput =
    """
      |task hello {
      |  File inputFile
      |
      |  command {
      |    echo ${inputFile}
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |  RUNTIME
      |}
      |
      |workflow hello {
      |  call hello
      |}
    """.stripMargin

  private val backendConfig = ConfigFactory.parseString(
    s"""{
        |  root = "local-cromwell-executions"
        |  filesystems {
        |    local {
        |      localization = [
        |        "hard-link", "soft-link", "copy"
        |      ]
        |    }
        |  }
        |}
        """.stripMargin)

  private val timeout = Timeout(1.seconds)

  after {
    Mockito.reset(sparkProcess)
  }

  override def afterAll(): Unit = system.shutdown()

  "executeTask method" should {
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
        response shouldBe a[SucceededResponse]
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
        response shouldBe a[FailedNonRetryableResponse]
        assert(response.asInstanceOf[FailedNonRetryableResponse].throwable.getMessage.contains(s"Execution process failed. Spark returned non zero status code:"))
      }

      cleanUpJob(jobPaths)
    }

    "return a successful task status even with a non-zero process exit" in {
      val runtime =
        """
          |runtime {
          | continueOnReturnCode: [911]
          |}
        """.stripMargin
      val jobDescriptor = prepareJob(runtimeString = runtime)
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
      when(stubProcess.exitValue()).thenReturn(911)
      when(sparkProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkProcess.processStderr).thenReturn(stderrResult)

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[SucceededResponse]
      }

      cleanUpJob(jobPaths)
    }

    "return a successful task status when it runs a docker command with working and output directory" in {
      val runtime =
        """
          |runtime {
          | docker: "ubuntu/latest"
          | dockerWorkingDir: "/workingDir"
          | dockerOutputDir: "/outputDir"
          |}
        """.stripMargin
      val jsonInputFile = createCannedFile("testFile", "some content").toPath.toAbsolutePath.toString
      val inputs = Map(
        "inputFile" -> WdlFile(jsonInputFile)
      )
      val jobDescriptor = prepareJob(helloWorldWdlWithFileInput, runtime, Option(inputs))
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
      when(stubProcess.exitValue()).thenReturn(0)
      when(sparkProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(sparkProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(sparkProcess.processStderr).thenReturn(stderrResult)

      whenReady(backend.execute) { response =>
        response shouldBe a[SucceededResponse]
      }

      val bashScript = Source.fromFile(jobPaths.script.toFile).getLines.mkString

      assert(bashScript.contains("docker run -w /workingDir -v"))
      assert(bashScript.contains(s"${jobPaths.script.getParent.toString}:"))
      assert(bashScript.contains("/call-hello:/outputDir --rm ubuntu/latest echo"))

      cleanUpJob(jobPaths)
    }
  }

  private def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): File = {
    val suffix = ".out"
    val file = dir match {
      case Some(path) => Files.createTempFile(path, prefix, suffix)
      case None => Files.createTempFile(prefix, suffix)
    }
    write(file.toFile, contents)
  }

  private def write(file: File, contents: String) = {
    val writer = new FileWriter(file)
    writer.write(contents)
    writer.flush()
    writer.close()
    file
  }

  private def cleanUpJob(jobPaths: JobPaths): Unit = jobPaths.workflowRoot.delete(true)

  private def prepareJob(source: String = helloWorldWdl, runtimeString: String = "", inputFiles: Option[Map[String, WdlValue]] = None): TestJobDescriptor = {
    val backendWorkflowDescriptor = buildWorkflowDescriptor(wdl = source, inputs = inputFiles.getOrElse(Map.empty), runtime = runtimeString)
    val backendConfigurationDescriptor = BackendConfigurationDescriptor(backendConfig, ConfigFactory.load)
    val jobDesc = jobDescriptorFromSingleCallWorkflow(backendWorkflowDescriptor, inputFiles.getOrElse(Map.empty))
    val jobPaths = new JobPaths(backendWorkflowDescriptor, backendConfig, jobDesc.key)
    val executionDir = jobPaths.callRoot
    val stdout = Paths.get(executionDir.path.toString, "stdout")
    stdout.toString.toFile.createIfNotExists(false)
    val submitFileStderr = executionDir.resolve("process.stderr")
    val submitFileStdout = executionDir.resolve("process.stdout")
    submitFileStdout.toString.toFile.createIfNotExists(false)
    submitFileStderr.toString.toFile.createIfNotExists(false)
    TestJobDescriptor(jobDesc, jobPaths, backendConfigurationDescriptor)
  }

  private case class TestJobDescriptor(jobDescriptor: BackendJobDescriptor, jobPaths: JobPaths, backendConfigurationDescriptor: BackendConfigurationDescriptor)

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
