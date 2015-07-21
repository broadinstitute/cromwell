package cromwell

import java.io.ByteArrayOutputStream

import akka.actor.ActorSystem
import akka.testkit.EventFilter
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.local.LocalBackend
import cromwell.server.WorkflowManagerSystem
import cromwell.util.SampleWdl.ThreeStep
import cromwell.util.{FileUtil, SampleWdl}
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.Await
import scala.concurrent.duration.Duration

class TestWorkflowManagerSystem extends WorkflowManagerSystem {
  override lazy val backend = new LocalBackend
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
      Main.highlight(Array(wdl))
    }
    val expected =
      s"""\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mps\u001b[0m {
         |  command {
         |    ps
         |  }
         |  output {
         |    \u001b[38;5;33mFile\u001b[0m \u001b[38;5;112mprocs\u001b[0m = \033[38;5;13mstdout\033[0m()
         |  }
         |}
         |
         |\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mcgrep\u001b[0m {
         |  command {
         |    grep '$${\u001b[38;5;33mString\u001b[0m \u001b[38;5;112mpattern\u001b[0m}' $${\u001b[38;5;33mFile\u001b[0m \u001b[38;5;112min_file\u001b[0m} | wc -l
         |  }
         |  output {
         |    \u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mcount\u001b[0m = \u001b[38;5;13mread_int\u001b[0m(\033[38;5;13mstdout\033[0m())
         |  }
         |}
         |
         |\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mwc\u001b[0m {
         |  command {
         |    cat $${\u001b[38;5;33mFile\u001b[0m \u001b[38;5;112min_file\u001b[0m} | wc -l
         |  }
         |  output {
         |    \u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mcount\u001b[0m = \u001b[38;5;13mread_int\u001b[0m(\033[38;5;13mstdout\033[0m())
         |  }
         |}
         |
         |\u001b[38;5;214mworkflow\u001b[0m \u001b[38;5;253mthree_step\u001b[0m {
         |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mps\u001b[0m
         |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mcgrep\u001b[0m {
         |    input: in_file=ps.procs
         |  }
         |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mwc\u001b[0m {
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
