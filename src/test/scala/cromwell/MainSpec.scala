package cromwell

import java.io.File

import cromwell.binding.WdlBinding
import cromwell.util.SampleWdl.ThreeStep
import cromwell.util.{FileUtil, SampleWdl}
import org.scalatest.{FlatSpec, Matchers}

class MainSpec extends FlatSpec with Matchers {
  val wdlFilePathAndWriter = FileUtil.tempFileAndWriter("wdl")
  val inputsJsonPathAndWriter = FileUtil.tempFileAndWriter("inputs")

  wdlFilePathAndWriter match {
    case (path, writer) =>
      writer.write(ThreeStep.WdlSource)
      writer.close()
  }

  inputsJsonPathAndWriter match {
    case (path, writer) =>
      writer.write("""{"three_step.cgrep.pattern": "..."}""")
      writer.close()
  }

  "Main" should "validate" in {
    val stream = new java.io.ByteArrayOutputStream()
    Console.withOut(stream) {
      Main.validate(Array(wdlFilePathAndWriter._1.toFile.getAbsolutePath))
    }
    stream.toString.length shouldEqual 0
  }

  it should "run" in {
    val stream = new java.io.ByteArrayOutputStream()
    Console.withOut(stream) {
      Main.run(Array(wdlFilePathAndWriter._1.toAbsolutePath.toString, inputsJsonPathAndWriter._1.toAbsolutePath.toString))
    }
    assert(stream.toString.contains("\"three_step.ps.procs\""))
    assert(stream.toString.contains("\"three_step.cgrep.count\""))
    assert(stream.toString.contains("\"three_step.wc.count\""))
  }

  it should "parse" in {
    val stream = new java.io.ByteArrayOutputStream()
    Console.withOut(stream) {
      Main.parse(Array(wdlFilePathAndWriter._1.toAbsolutePath.toString))
    }
    assert(stream.toString.contains("(Document:"))
  }

  it should "return inputs" in {
    val stream = new java.io.ByteArrayOutputStream()
    Console.withOut(stream) {
      Main.inputs(Array(wdlFilePathAndWriter._1.toAbsolutePath.toString))
    }
    assert(stream.toString.contains("\"three_step.cgrep.pattern\""))
  }

  it should "print usage" in {
    val stream = new java.io.ByteArrayOutputStream()
    Console.withOut(stream) {
      Main.usageAndExit(false)
    }
    assert(stream.toString.contains("java -jar cromwell.jar <action> <parameters>"))
  }
}
