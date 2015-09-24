package cromwell

import java.io.ByteArrayOutputStream

import akka.actor.ActorSystem
import akka.testkit.EventFilter
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.local.LocalBackend
import cromwell.parser.BackendType
import cromwell.server.WorkflowManagerSystem
import cromwell.util.SampleWdl.ThreeStep
import cromwell.util.{FileUtil, SampleWdl}
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.Await
import scala.concurrent.duration.Duration

class TestWorkflowManagerSystem extends WorkflowManagerSystem {
  override lazy val backendType = BackendType.LOCAL
  override implicit val actorSystem = ActorSystem(systemName, ConfigFactory.parseString(CromwellTestkitSpec.ConfigText))
}

class MainSpec extends FlatSpec with Matchers {

  private def wdlAndInputs(sampleWdl: SampleWdl): (String, String) = {
    val wdlFilePathAndWriter = FileUtil.tempFileAndWriter("wdl")
    val inputsJsonPathAndWriter = FileUtil.tempFileAndWriter("inputs")

    wdlFilePathAndWriter match {
      case (path, writer) =>
        writer.write(sampleWdl.wdlSource())
        writer.close()
    }
    inputsJsonPathAndWriter match {
      case (path, writer) =>
        writer.write(sampleWdl.wdlJson)
        writer.close()
    }

    (wdlFilePathAndWriter._1.toFile.getAbsolutePath, inputsJsonPathAndWriter._1.toFile.getAbsolutePath)
  }

  def baos = new ByteArrayOutputStream()

  "Main" should "validate" in {
    val stream = baos
    Console.withOut(stream) {
      val (wdl, _) = wdlAndInputs(ThreeStep)
      Main.validate(Array(wdl))
    }
    stream.toString.length shouldEqual 0
  }

  it should "parse" in {
    val stream = baos
    Console.withOut(stream) {
      val (wdl, _) = wdlAndInputs(ThreeStep)
      Main.parse(Array(wdl))
    }
    assert(stream.toString.contains("(Document:"))
  }

  it should "highlight" in {
    val stream = baos
    Console.withOut(stream) {
      val (wdl, _) = wdlAndInputs(ThreeStep)
      Main.highlight(Array(wdl, "html"))
    }
    val expected =
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
    stream.toString.stripLineEnd shouldEqual expected
  }

  it should "return inputs" in {
    val stream = baos
    Console.withOut(stream) {
      val (wdl, _) = wdlAndInputs(ThreeStep)
      Main.inputs(Array(wdl))
    }
    assert(stream.toString.contains("\"three_step.cgrep.pattern\""))
  }

  it should "run" in {
    val workflowManagerSystem = new TestWorkflowManagerSystem
    implicit val system = workflowManagerSystem.actorSystem
    EventFilter.info(pattern = s"transitioning from Running to Succeeded.", occurrences = 1).intercept {
      val (wdl, inputs) = wdlAndInputs(ThreeStep)
      Main.run(Array(wdl, inputs), workflowManagerSystem)
    }
    Await.result(workflowManagerSystem.shutdown(), Duration.Inf)
  }

  it should "print usage" in {
    val stream = baos
    Console.withOut(stream) {
      Main.usageAndExit(false)
    }
    assert(stream.toString.contains("java -jar cromwell.jar <action> <parameters>"))
  }
}
