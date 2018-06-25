package womtool

import better.files._
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import womtool.WomtoolMainSpec._
import womtool.SampleWdl.{EmptyTask, EmptyWorkflow, ThreeStep}

class WomtoolMainSpec extends FlatSpec with Matchers with BeforeAndAfterAll {

  import WomtoolMain._

  behavior of "Main"

  val threeStep = ThreeStep.wdlSource()

  it should "print usage" in {
    WomtoolMain.runWomtool(Seq.empty[String]) match{
      case BadUsageTermination(msg) => msg should include("Usage: java -jar womtool.jar [validate|inputs|parse|highlight|graph|womgraph] [options] workflow-source")
      case other => fail(s"Expected BadUsageTermination but got $other")
    }
  }

  it should "parse" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      val res = WomtoolMain.runWomtool(Seq("parse", wdlAndInputs.wdl))
      assert(res.isInstanceOf[SuccessfulTermination])
      res.stdout.get should include("(Namespace:")
    }
  }

  it should "highlight" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      val res = WomtoolMain.runWomtool(Seq("highlight", "-h", "html", wdlAndInputs.wdl))
      assert(res.isInstanceOf[SuccessfulTermination])
      res.stdout.get.stripLineEnd should be(HighlightedWdlHtml)
    }
  }

  it should "highlight using console highlighting" in {
    testWdl(EmptyWorkflow) { wdlAndInputs =>
      val res = WomtoolMain.runWomtool(Seq("highlight", wdlAndInputs.wdl, "--highlight-mode", "console"))
      assert(res.isInstanceOf[SuccessfulTermination])
      res.stdout.get.stripLineEnd should include("empty_workflow")
    }
  }

  it should "return inputs" in {
    testWdl(ThreeStep) { wdlAndInputs =>
      val res = WomtoolMain.runWomtool(Seq("inputs", wdlAndInputs.wdl))
      assert(res.isInstanceOf[SuccessfulTermination])
      res.stdout.get should include("\"three_step.cgrep.pattern\"")
    }
  }

  it should "not return inputs when there is no workflow" in {
    testWdl(EmptyTask) { wdlAndInputs =>
      val res = WomtoolMain.runWomtool(Seq("inputs", wdlAndInputs.wdl))
      res should be(UnsuccessfulTermination("Cannot convert WOM bundle to executable. No primary callable was available."))
    }
  }
}

object WomtoolMainSpec {
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

  val HighlightedWdlHtml =
    """<span class="keyword">task</span> <span class="name">ps</span> {
      |  <span class="section">command</span> {
      |    <span class="command">ps</span>
      |  }
      |  <span class="section">output</span> {
      |    <span class="type">File</span> <span class="variable">procs</span> = <span class="function">stdout</span>()
      |  }
      |}
      |
      |<span class="keyword">task</span> <span class="name">cgrep</span> {
      |  <span class="type">String</span> <span class="variable">pattern</span>
      |  <span class="type">File</span> <span class="variable">in_file</span>
      |  <span class="section">command</span> {
      |    <span class="command">grep '${pattern}' ${in_file} | wc -l</span>
      |  }
      |  <span class="section">output</span> {
      |    <span class="type">Int</span> <span class="variable">count</span> = <span class="function">read_int</span>(<span class="function">stdout</span>())
      |  }
      |}
      |
      |<span class="keyword">task</span> <span class="name">wc</span> {
      |  <span class="type">File</span> <span class="variable">in_file</span>
      |  <span class="section">command</span> {
      |    <span class="command">cat ${in_file} | wc -l</span>
      |  }
      |  <span class="section">output</span> {
      |    <span class="type">Int</span> <span class="variable">count</span> = <span class="function">read_int</span>(<span class="function">stdout</span>())
      |  }
      |}
      |
      |<span class="keyword">workflow</span> <span class="name">three_step</span> {
      |  <span class="keyword">call</span> <span class="name">ps</span>
      |  <span class="keyword">call</span> <span class="name">cgrep</span> {
      |    input: in_file=ps.procs
      |  }
      |  <span class="keyword">call</span> <span class="name">wc</span> {
      |    input: in_file=ps.procs
      |  }
      |}""".stripMargin

  /**
    * Create a temporary wdl file and inputs for the sampleWdl.
    * When the various properties are lazily accessed, they are also registered for deletion after the suite completes.
    */
  case class WdlAndInputs(sampleWdl: SampleWdl, optionsJson: String = "{}") {
    // Track all the temporary files we create, and delete them after the test.
    private var tempFiles = Vector.empty[File]

    lazy val wdlFile = {
      val path = File.newTemporaryFile(s"${sampleWdl.name}.", ".wdl")
      tempFiles :+= path
      path write sampleWdl.wdlSource("")
      path
    }

    lazy val wdl = wdlFile.pathAsString

    lazy val inputsFile = {
      val path = swapExt(wdlFile, ".wdl", ".inputs")
      tempFiles :+= path
      path write sampleWdl.wdlJson
      path
    }

    lazy val inputs = inputsFile.pathAsString

    lazy val optionsFile = {
      val path = swapExt(wdlFile, ".wdl", ".options")
      tempFiles :+= path
      path write optionsJson
      path
    }

    lazy val options = optionsFile.pathAsString

    lazy val metadataFile = {
      val path = swapExt(wdlFile, ".wdl", ".metadata.json")
      tempFiles :+= path
      path
    }

    lazy val metadata = metadataFile.pathAsString

    def deleteTempFiles() = tempFiles.foreach(_.delete(swallowIOExceptions = true))
  }

  def swapExt(filePath: File, oldExt: String, newExt: String): File = {
    File(filePath.toString.stripSuffix(oldExt) + newExt)
  }
}
