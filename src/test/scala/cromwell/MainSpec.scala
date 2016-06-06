package cromwell

import java.io.{ByteArrayOutputStream, OutputStream}
import java.nio.file.Path
import java.text.SimpleDateFormat
import java.util.Date

import akka.util.Timeout
import better.files._
import cromwell.util.FileUtil._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl._
import org.apache.commons.io.output.TeeOutputStream
import org.scalatest.concurrent.TimeLimitedTests
import org.scalatest.time.Span
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class MainSpec extends FlatSpec with Matchers with BeforeAndAfterAll with TimeLimitedTests {
  import MainSpec._

  behavior of "Main"

  override val timeLimit: Span = 60 seconds

  it should "print usage" in {
    assert(traceMain(_.usageAndExit()).out.contains(UsageSnippet))
  }

  it should "print usage when no args" in {
    val result = traceAction()
    assert(result.out.contains(UsageSnippet))
  }

  it should "run" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      val wdl = wdlAndInputs.wdl
      val inputs = wdlAndInputs.inputs
      traceInfoRun(wdl, inputs)("transitioning from Running to Succeeded.") should be(0)
    }
  }

  it should "run using args" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      val wdl = wdlAndInputs.wdl
      val inputs = wdlAndInputs.inputs
      traceInfoAction("run", wdl, inputs)("transitioning from Running to Succeeded.") should be(0)
    }
  }

  it should "run and locate the default inputs path" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      val wdl = wdlAndInputs.wdl
      wdlAndInputs.inputs
      traceInfoRun(wdl)("transitioning from Running to Succeed.") should be(0)
    }
  }

  it should "run if the inputs path is \"-\"" in {
    testWdl(GoodbyeWorld) { wdlAndInputs =>
      val wdl = wdlAndInputs.wdl
      traceErrorWithExceptionRun(wdl, "-", "-", "-")("transitioned to state Failed") should be(1)
    }
  }

  it should "run reading options" in {
    testWdl(ThreeStep, optionsJson = """{ foobar bad json! }""") { wdlAndInputs =>
      val wdl = wdlAndInputs.wdl
      val inputs = wdlAndInputs.inputs
      val options = wdlAndInputs.options
      traceErrorWithExceptionRun(wdl, inputs, options)(".*Workflow failed submission:", classOf[IllegalArgumentException]) should be(1)
    }
  }

  it should "run writing metadata" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      val wdl = wdlAndInputs.wdl
      val inputs = wdlAndInputs.inputs
      val metadata = wdlAndInputs.metadata
      traceInfoRun(wdl, inputs, "-", metadata)("transitioning from Running to Succeeded.") should be(0)
      assert(wdlAndInputs.metadataPath.contentAsString.contains("\"three_step.cgrep.pattern\""))
    }
  }

  it should "fail run if the inputs path is not found, and not set to -" in {
    testWdl(EmptyWorkflow) { wdlAndInputs =>
      val wdl = wdlAndInputs.wdl
      val result = traceMain(_.run(Array(wdl)))
      assert(result.err.contains("Inputs does not exist"))
      result.returnCode should be(1)
    }
  }

  it should "fail run if the inputs path is not readable" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      wdlAndInputs.inputsPath setPermissions Set.empty
      val wdl = wdlAndInputs.wdl
      val inputs = wdlAndInputs.inputs
      val result = traceMain(_.run(Array(wdl, inputs)))
      assert(result.err.contains("Inputs is not readable"))
      result.returnCode should be(1)
    }
  }

  it should "fail run if the inputs path is not valid inputs json" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      wdlAndInputs.inputsPath write "[]"
      val wdl = wdlAndInputs.wdl
      val inputs = wdlAndInputs.inputs
      val result = traceMain(_.run(Array(wdl, inputs)))
      assert(result.err.contains("Workflow inputs JSON cannot be parsed to JsObject"))
      result.returnCode should be(1)
    }
  }

  it should "fail run with not enough args" in {
    val result = traceMain(_.run(Array.empty[String]))
    assert(result.out.contains(UsageSnippet))
    result.returnCode should be(-1)
  }

  it should "fail run with too many args" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      val wdl = wdlAndInputs.wdl
      val inputs = wdlAndInputs.inputs
      val result = traceMain(_.run(Array(wdl, inputs, "-", "-", "-", "extra")))
      assert(result.out.contains(UsageSnippet))
      result.returnCode should be(-1)
    }
  }
}

object MainSpec {
  import CromwellTestkitSpec._

  implicit val AskTimeout = Timeout(5 seconds)

  /** The return code, plus any captured text from Console.stdout and Console.stderr while executing a block. */
  case class TraceResult(returnCode: Int, out: String, err: String)
  private val dateFormat = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS")
  private def now = dateFormat.format(new Date)

  /**
    * Tests running a sample wdl, providing the inputs, and cleaning up the temp files only if no exceptions occur.
    *
    * @param sampleWdl The sample wdl to run.
    * @param optionsJson Optional json for the options file.
    * @param block The block provided the inputs, returning some value.
    * @tparam T The return type of the block.
    * @return The result of running the block.
    */
  def testWdl[T](sampleWdl: SampleWdl, optionsJson: String = "{}")(block: WdlAndInputs => T): T = {
    val wdlAndInputs = WdlAndInputs(sampleWdl, optionsJson)
    val result = block(wdlAndInputs)
    wdlAndInputs.deleteTempFiles()
    result
  }

  /**
    * Saves and restores some system properties.
    *
    * @param keys SystemProperty names that should be saved.
    * @param block The block to run.
    * @tparam T The return type of the block.
    * @return The result of running the block.
    */
  def modifyingSysProps[T](keys: String*)(block: => T): T = {
    val saved = sys.props.filterKeys(keys.contains).toMap
    try {
      block
    } finally {
      val restore = sys.props
      keys.foreach(restore.-=)
      saved.foreach(restore.+=)
    }
  }

  /**
    * Prints the entry and exit of a block.
    *
    * Used primarily to (hopefully?) force a heisenbug to hide itself.
    *
    * The pattern matching utilities wait for a certain event to occur within a block. All system events are also
    * println'ed to the output. Each block supposedly also runs system.awaitTermination() at the end of Main.run,
    * blocking until the entire system is shutdown.
    *
    * However, on Travis, a number of failures have been seen where the pattern appears in the console output, but the
    * filter returns false that the pattern had been seen.
    *
    * println is internally synchronized, so it's possible adding a call to printBlock _inside_ the waitFor... is masking
    * a heisenbug, allowing the other events to be processed semi-synchronously before the final println in printBlock is
    * allowed to finish. Switching the order of waitFor... and printBlock _seemed_ to reduce the intermittent errors,
    * but this hasn't been exhaustively confirmed. Looking forward to seeing the WorkflowActor system revamp its concept
    * of "terminated" at some point in the future anyway.
    */
  private def printBlock(action: String, args: Seq[String])(block: => Int): Int = {
    try {
      println(s"[$now] [block] Entering block: $action ${args.mkString(" ")}")
      block
    } finally {
      println(s"[$now] [block] Exiting block: $action ${args.mkString(" ")}")
    }
  }

  /**
    * Runs the "run" method and waits for a particular pattern to appear as a Log Info event in the system.
    *
    * @param args Args to pass to Main.run().
    * @param pattern The pattern to watch for.
    * @return The return code of run.
    */
  def traceInfoRun(args: String*)(pattern: String): Int = {
    withTestWorkflowManagerSystem { workflowManagerSystem =>
      waitForInfo(pattern)(
        printBlock("run", args) {
          new Main(workflowManagerSystem).run(args)
        }
      )(workflowManagerSystem.actorSystem)
    }
  }

  /**
    * Runs the "run" method and waits for a particular pattern to appear as a Log Error event in the system, with an
    * attached exception.
    *
    * @param args Args to pass to Main.run().
    * @param pattern The pattern to watch for.
    * @return The return code of run.
    */
  def traceErrorWithExceptionRun(args: String*)(pattern: String, throwableClass: Class[_ <: Throwable] = classOf[Throwable]): Int = {
    withTestWorkflowManagerSystem { workflowManagerSystem =>
      waitForErrorWithException(pattern, throwableClass = throwableClass)(
        printBlock("run", args) {
          new Main(workflowManagerSystem).run(args)
        }
      )(workflowManagerSystem.actorSystem)
    }
  }

  /**
    * Runs the "runAction" method and waits for a particular pattern to appear as a Log Info event in the system.
    *
    * @param args Args to pass to Main.run().
    * @param pattern The pattern to watch for
    * @return The return code of run.
    */
  def traceInfoAction(args: String*)(pattern: String): Int = {
    withTestWorkflowManagerSystem { workflowManagerSystem =>
      waitForInfo(pattern)(
        printBlock("runAction", args) {
          new Main(workflowManagerSystem).runAction(args) match {
            case status: Int => status
          }
        }
      )(workflowManagerSystem.actorSystem)
    }
  }

  /**
    * Loans an instance of Main, returning the return code plus everything that passed through Console.out/Console.err
    * during the run.
    *
    * @param block Block to run.
    * @return return code plus Console.out/Console.err during the block.
    */
  def traceMain(block: Main => Int): TraceResult = {
    withTestWorkflowManagerSystem { workflowManagerSystem =>
      val outStream = TeeStream(Console.out)
      val errStream = TeeStream(Console.err)
      val status = {
        Console.withOut(outStream.teed) {
          Console.withErr(errStream.teed) {
            block(new Main(workflowManagerSystem))
          }
        }
      }
      TraceResult(status, outStream.captured, errStream.captured)
    }
  }

  /**
    * Runs Main.runAction, returning the return code plus everything that passed through Console.out/Console.err
    * during the run.
    *
    * @param args Arguments to pass to runAction.
    * @return return code plus Console.out/Console.err during the block.
    */
  def traceAction(args: String*): TraceResult = {
    traceMain { main =>
      main.runAction(args) match {
        case status: Int => status
      }
    }
  }

  /**
    * Create a temporary wdl file and inputs for the sampleWdl.
    * When the various properties are lazily accessed, they are also registered for deletion after the suite completes.
    */
  case class WdlAndInputs(sampleWdl: SampleWdl, optionsJson: String = "{}") {
    // Track all the temporary files we create, and delete them after the test.
    private var tempFiles = Vector.empty[Path]

    lazy val wdlPath: Path = {
      val path = File.newTemp(s"${sampleWdl.name}.", ".wdl").path
      tempFiles :+= path
      path write sampleWdl.wdlSource("")
      path
    }

    lazy val wdl = wdlPath.fullPath

    lazy val inputsPath = {
      val path = wdlPath.swapExt(".wdl", ".inputs")
      tempFiles :+= path
      path write sampleWdl.wdlJson
      path
    }

    lazy val inputs = inputsPath.fullPath

    lazy val optionsPath = {
      val path = wdlPath.swapExt(".wdl", ".options")
      tempFiles :+= path
      path write optionsJson
      path
    }

    lazy val options = optionsPath.fullPath

    lazy val metadataPath = {
      val path = wdlPath.swapExt(".wdl", ".metadata.json")
      tempFiles :+= path
      path.toAbsolutePath
    }

    lazy val metadata = metadataPath.fullPath

    def deleteTempFiles() = tempFiles.foreach(_.delete(ignoreIOExceptions = true))
  }

  /**
    * Utility for capturing output while also streaming it to stdout/stderr.
    *
    * @param orig The stream to share.
    */
  case class TeeStream(orig: OutputStream) {
    private lazy val byteStream = new ByteArrayOutputStream()

    /** The teed stream. One should NOT call close on the stream, as it combo of a byte stream and a console stream. */
    lazy val teed = new TeeOutputStream(orig, byteStream)

    /** The captured text from the original stream. */
    def captured = byteStream.toString
  }

  val UsageSnippet = "java -jar cromwell.jar <action> <parameters>"
}
