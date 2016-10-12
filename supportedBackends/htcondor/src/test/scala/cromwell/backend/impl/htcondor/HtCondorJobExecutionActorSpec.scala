package cromwell.backend.impl.htcondor

import java.io.Writer
import java.nio.file.{Files, Path}

import akka.actor.{Actor, Props}
import akka.testkit.{ImplicitSender, TestActorRef}
import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.{JobFailedNonRetryableResponse, JobSucceededResponse}
import cromwell.backend.impl.htcondor.caching.CacheActor
import cromwell.backend.impl.htcondor.caching.exception.CachedResultNotFoundException
import cromwell.backend.impl.htcondor.caching.model.CachedExecutionResult
import cromwell.backend.io.JobPathsWithDocker
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendSpec}
import cromwell.core._
import cromwell.core.path.{PathWriter, TailedWriter, UntailedWriter}
import cromwell.services.keyvalue.KeyValueServiceActor.{KvGet, KvPair, KvPut}
import org.mockito.Mockito
import org.mockito.Mockito._
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.mockito.MockitoSugar
import org.scalatest.{BeforeAndAfter, Matchers, WordSpecLike}
import wdl4s.types.{WdlArrayType, WdlFileType}
import wdl4s.values.{WdlArray, WdlFile, WdlValue}

import scala.concurrent.duration._
import scala.io.Source
import scala.sys.process.{Process, ProcessLogger}

class HtCondorJobExecutionActorSpec extends TestKitSuite("HtCondorJobExecutionActorSpec")
  with WordSpecLike
  with Matchers
  with MockitoSugar
  with BeforeAndAfter
  with ImplicitSender {

  import BackendSpec._

  private val htCondorCommands: HtCondorCommands = new HtCondorCommands
  private val htCondorProcess: HtCondorProcess = mock[HtCondorProcess]
  private val cacheActorMockProps = Props(new CacheActorMock())

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
      |workflow wf_hello {
      |  call hello
      |}
    """.stripMargin

  private val helloWorldWdlWithFileInput =
    s"""
      |task hello {
      |  File inputFile
      |
      |  command {
      |    echo $${inputFile}
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

  private val helloWorldWdlWithFileArrayInput =
    s"""
      |task hello {
      |  Array[File] inputFiles
      |
      |  command {
      |    echo $${sep=' ' inputFiles}
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

  private val backendConfig = ConfigFactory.parseString(
    s"""{
        |  root = "local-cromwell-executions"
        |
        |  docker {
        |    cmd = "docker run -w %s %s %s %s --rm %s /bin/bash -c \\"%s\\""
        |    defaultWorkingDir = "/workingDir/"
        |    defaultOutputDir = "/output/"
        |  }
        |
        |  filesystems {
        |    local {
        |      localization = [
        |        "hard-link", "soft-link", "copy"
        |      ]
        |    }
        |  }
        |  poll-interval = 3
        |}
        """.stripMargin)

  private val timeout = Timeout(1.seconds)

  after {
    Mockito.reset(htCondorProcess)
  }

  "executeTask method" should {
    "return succeeded task status with stdout" in {
      val jobDescriptor = prepareJob()
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""

      when(htCondorProcess.commandList(any[String])).thenReturn(Seq.empty[String])
      when(htCondorProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(htCondorProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(htCondorProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(htCondorProcess.processStderr).thenReturn(stderrResult)
      when(htCondorProcess.jobReturnCode(any[String], any[Path])).thenReturn(Option(0))

      val backend = TestActorRef(new HtCondorJobExecutionActor(job, backendConfigDesc, system.deadLetters, None) {
        override lazy val cmds = htCondorCommands
        override lazy val extProcess = htCondorProcess
      }).underlyingActor

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobSucceededResponse]
        verify(htCondorProcess, times(1)).externalProcess(any[Seq[String]], any[ProcessLogger])
        verify(htCondorProcess, times(1)).tailedWriter(any[Int], any[Path])
        verify(htCondorProcess, times(1)).untailedWriter(any[Path])
      }

      cleanUpJob(jobPaths)
    }

    "return succeeded task status when it recovers from a shutdown" in {
      val jobDescriptor = prepareJob()
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""
      val kVServiceActor = system.actorOf(Props(new KVServiceActor()))

      when(htCondorProcess.commandList(any[String])).thenReturn(Seq.empty[String])
      when(htCondorProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(htCondorProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(htCondorProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(htCondorProcess.processStderr).thenReturn(stderrResult)
      when(htCondorProcess.jobReturnCode(any[String], any[Path])).thenReturn(Option(0))

      val backend = TestActorRef(new HtCondorJobExecutionActor(job, backendConfigDesc, kVServiceActor, None) {
        override lazy val cmds = htCondorCommands
        override lazy val extProcess = htCondorProcess
      }).underlyingActor

      whenReady(backend.recover, timeout) { response =>
        response shouldBe a[JobSucceededResponse]
      }

      cleanUpJob(jobPaths)
    }

    "return succeeded task status with stdout when cache is enabled" in {
      val jobDescriptor = prepareJob()
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""

      when(htCondorProcess.commandList(any[String])).thenReturn(Seq.empty[String])
      when(htCondorProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(htCondorProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(htCondorProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(htCondorProcess.processStderr).thenReturn(stderrResult)
      when(htCondorProcess.jobReturnCode(any[String], any[Path])).thenReturn(Option(0))

      val backend = TestActorRef(new HtCondorJobExecutionActor(job, backendConfigDesc, system.deadLetters, Some(cacheActorMockProps)) {
        override lazy val cmds = htCondorCommands
        override lazy val extProcess = htCondorProcess
      }).underlyingActor

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobSucceededResponse]
        verify(htCondorProcess, times(1)).externalProcess(any[Seq[String]], any[ProcessLogger])
        verify(htCondorProcess, times(1)).tailedWriter(any[Int], any[Path])
        verify(htCondorProcess, times(1)).untailedWriter(any[Path])
      }

      cleanUpJob(jobPaths)
    }

    "return failed task status with stderr on non-zero process exit" in {
      val jobDescriptor = prepareJob()
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)

      val backend = TestActorRef(new HtCondorJobExecutionActor(job, backendConfigDesc, system.deadLetters, Some(cacheActorMockProps)) {
        override lazy val cmds = htCondorCommands
        override lazy val extProcess = htCondorProcess
      }).underlyingActor
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""

      when(htCondorProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(htCondorProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(htCondorProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(htCondorProcess.processStderr).thenReturn(stderrResult)
      when(htCondorProcess.jobReturnCode(any[String], any[Path])).thenReturn(Option(-1))

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains("Job exited with invalid return code"))
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

      val backend = TestActorRef(new HtCondorJobExecutionActor(job, backendConfigDesc, system.deadLetters, Some(cacheActorMockProps)) {
        override lazy val cmds = htCondorCommands
        override lazy val extProcess = htCondorProcess
      }).underlyingActor
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""

      when(htCondorProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(htCondorProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(htCondorProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(htCondorProcess.processStderr).thenReturn(stderrResult)
      when(htCondorProcess.jobReturnCode(any[String], any[Path])).thenReturn(Option(911))

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobSucceededResponse]
      }

      cleanUpJob(jobPaths)
    }

    "return a successful task status when it runs a docker command with working and output directory" in {
      val runtime =
        """
          |runtime {
          | docker: "ubuntu/latest"
          | dockerWorkingDir: "/workingDir/"
          | dockerOutputDir: "/outputDir/"
          |}
        """.stripMargin
      val jsonInputFile = createCannedFile("testFile", "some content").pathAsString
      val inputs = Map(
        "wf_hello.hello.inputFile" -> WdlFile(jsonInputFile)
      )
      val jobDescriptor = prepareJob(helloWorldWdlWithFileInput, runtime, Option(inputs))
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)

      val backend = TestActorRef(new HtCondorJobExecutionActor(job, backendConfigDesc, system.deadLetters, Some(cacheActorMockProps)) {
        override lazy val cmds = htCondorCommands
        override lazy val extProcess = htCondorProcess
      }).underlyingActor
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""

      when(htCondorProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(htCondorProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(htCondorProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(htCondorProcess.processStderr).thenReturn(stderrResult)
      when(htCondorProcess.jobReturnCode(any[String], any[Path])).thenReturn(Option(0))

      whenReady(backend.execute) { response =>
        response shouldBe a[JobSucceededResponse]
      }

      val bashScript = Source.fromFile(jobPaths.script.toFile).getLines.mkString

      assert(bashScript.contains("docker run -w /workingDir/ -v"))
      assert(bashScript.contains(":/workingDir/"))
      assert(bashScript.contains(":ro"))
      assert(bashScript.contains("/call-hello/execution:/outputDir/ --rm ubuntu/latest /bin/bash -c \"echo"))

      cleanUpJob(jobPaths)
    }

    "return failed when cmds fails to write script" in {
      val htCondorCommandsMock: HtCondorCommands = mock[HtCondorCommands]
      val jobDescriptor = prepareJob()
      val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)

      val backend = TestActorRef(new HtCondorJobExecutionActor(job, backendConfigDesc, system.deadLetters, Some(cacheActorMockProps)) {
        override lazy val cmds = htCondorCommandsMock
        override lazy val extProcess = htCondorProcess
      }).underlyingActor
      val stubProcess = mock[Process]
      val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
      val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
      val stderrResult = ""

      when(htCondorCommandsMock.writeScript(any[String], any[Path], any[Path])).thenThrow(new IllegalStateException("Could not write the file."))
      when(htCondorProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
      when(stubProcess.exitValue()).thenReturn(0)
      when(htCondorProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
      when(htCondorProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
      when(htCondorProcess.processStderr).thenReturn(stderrResult)
      when(htCondorProcess.jobReturnCode(any[String], any[Path])).thenReturn(Option(-1))

      whenReady(backend.execute, timeout) { response =>
        response shouldBe a[JobFailedNonRetryableResponse]
        assert(response.asInstanceOf[JobFailedNonRetryableResponse].throwable.getMessage.contains("Could not write the file."))
      }

      cleanUpJob(jobPaths)
    }
  }

  "return a successful task status when it tries to run a docker command containing file data from a WDL file array" in {
    val runtime =
      """
        |runtime {
        | docker: "ubuntu/latest"
        | dockerWorkingDir: "/workingDir/"
        | dockerOutputDir: "/outputDir/"
        |}
      """.stripMargin

    val tempDir1 = Files.createTempDirectory("dir1")
    val tempDir2 = Files.createTempDirectory("dir2")
    val jsonInputFile =
      createCannedFile(prefix = "testFile", contents = "some content", dir = Some(tempDir1)).pathAsString
    val jsonInputFile2 =
      createCannedFile(prefix = "testFile2", contents = "some other content", dir = Some(tempDir2)).pathAsString

    val inputs = Map(
      "wf_hello.hello.inputFiles" -> WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile(jsonInputFile), WdlFile(jsonInputFile2)))
    )
    val jobDescriptor = prepareJob(helloWorldWdlWithFileArrayInput, runtime, Option(inputs))
    val (job, jobPaths, backendConfigDesc) = (jobDescriptor.jobDescriptor, jobDescriptor.jobPaths, jobDescriptor.backendConfigurationDescriptor)

    val backend = TestActorRef(new HtCondorJobExecutionActor(job, backendConfigDesc, system.deadLetters, Some(cacheActorMockProps)) {
      override lazy val cmds = htCondorCommands
      override lazy val extProcess = htCondorProcess
    }).underlyingActor
    val stubProcess = mock[Process]
    val stubUntailed = new UntailedWriter(jobPaths.stdout) with MockPathWriter
    val stubTailed = new TailedWriter(jobPaths.stderr, 100) with MockPathWriter
    val stderrResult = ""

    when(htCondorProcess.externalProcess(any[Seq[String]], any[ProcessLogger])).thenReturn(stubProcess)
    when(stubProcess.exitValue()).thenReturn(0)
    when(htCondorProcess.tailedWriter(any[Int], any[Path])).thenReturn(stubTailed)
    when(htCondorProcess.untailedWriter(any[Path])).thenReturn(stubUntailed)
    when(htCondorProcess.processStderr).thenReturn(stderrResult)
    when(htCondorProcess.jobReturnCode(any[String], any[Path])).thenReturn(Option(0))

    whenReady(backend.execute) { response =>
      response shouldBe a[JobSucceededResponse]
    }

    val bashScript = Source.fromFile(jobPaths.script.toFile).getLines.mkString

    assert(bashScript.contains("docker run -w /workingDir/ -v"))
    assert(bashScript.contains(":/workingDir/"))
    assert(bashScript.contains(tempDir1.toAbsolutePath.toString))
    assert(bashScript.contains(tempDir2.toAbsolutePath.toString))
    assert(bashScript.contains("/call-hello/execution:/outputDir/ --rm ubuntu/latest /bin/bash -c \"echo"))

    cleanUpJob(jobPaths)
  }

  private def cleanUpJob(jobPaths: JobPathsWithDocker): Unit = {
    File(jobPaths.workflowRoot).delete(true)
    ()
  }

  private def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): File = {
    val suffix = ".out"
    val file = File.newTemporaryFile(prefix, suffix, dir.map(File.apply))
    file.write(contents)
  }

  val emptyWorkflowOptions = WorkflowOptions.fromMap(Map.empty).get

  private def prepareJob(source: String = helloWorldWdl, runtimeString: String = "", inputFiles: Option[Map[String, WdlValue]] = None): TestJobDescriptor = {
    val backendWorkflowDescriptor = buildWorkflowDescriptor(wdl = source, inputs = inputFiles.getOrElse(Map.empty), runtime = runtimeString)
    val backendConfigurationDescriptor = BackendConfigurationDescriptor(backendConfig, ConfigFactory.load)
    val jobDesc = jobDescriptorFromSingleCallWorkflow(backendWorkflowDescriptor, inputFiles.getOrElse(Map.empty), emptyWorkflowOptions, Set.empty)
    val jobPaths = new JobPathsWithDocker(jobDesc.key, backendWorkflowDescriptor, backendConfig)
    val executionDir = File(jobPaths.callExecutionRoot)
    val stdout = File(executionDir.pathAsString, "stdout")
    stdout.createIfNotExists(asDirectory = false, createParents = true)
    val submitFileStderr = executionDir./("submitfile.stderr")
    val submitFileStdout = executionDir./("submitfile.stdout")
    submitFileStdout.createIfNotExists(asDirectory = false, createParents = true)
    submitFileStdout <<
      """Submitting job(s)..
        |1 job(s) submitted to cluster 88.
      """.stripMargin.trim
    submitFileStderr.createIfNotExists(asDirectory = false, createParents = true)
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

  class CacheActorMock extends CacheActor {
    override def readExecutionResult(hash: String): CachedExecutionResult = throw new CachedResultNotFoundException("Entry not found.")

    override def storeExecutionResult(cachedExecutionResult: CachedExecutionResult): Unit = ()
  }

  class KVServiceActor extends Actor {
    override def receive: Receive = {
      case KvPut => // Do nothing
      case KvGet(kvKey) => sender ! KvPair(kvKey, Option("123"))
    }
  }

}
