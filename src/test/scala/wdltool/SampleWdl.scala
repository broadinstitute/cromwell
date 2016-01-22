package wdltool

import java.io.{FileWriter, File}
import java.nio.file.{Files, Path}
import spray.json._

import wdl4s._
import wdl4s.values._

import scala.language.postfixOps

trait SampleWdl {
  def wdlSource(runtime: String = ""): WdlSource
  def rawInputs: WorkflowRawInputs

  def name = getClass.getSimpleName.stripSuffix("$")

  implicit object AnyJsonFormat extends JsonFormat[Any] {
    def write(x: Any) = x match {
      case n: Int => JsNumber(n)
      case s: String => JsString(s)
      case b: Boolean => if(b) JsTrue else JsFalse
      case s: Seq[Any] => JsArray(s map {_.toJson} toVector)
      case a: WdlArray => write(a.value)
      case s: WdlString => JsString(s.value)
      case i: WdlInteger => JsNumber(i.value)
      case f: WdlFloat => JsNumber(f.value)
      case f: WdlFile => JsString(f.value)
    }
    def read(value: JsValue) = throw new NotImplementedError(s"Reading JSON not implemented: $value")
  }

  implicit object RawInputsJsonFormat extends JsonFormat[WorkflowRawInputs] {
    def write(inputs: WorkflowRawInputs) = JsObject(inputs map { case (k, v) => k -> v.toJson })
    def read(value: JsValue) = throw new NotImplementedError(s"Reading JSON not implemented: $value")
  }

  def wdlJson: WdlJson = rawInputs.toJson.prettyPrint

  def createFileArray(base: Path): Unit = {
    createFile("f1", base, "line1\nline2\n")
    createFile("f2", base, "line3\nline4\n")
    createFile("f3", base, "line5\n")
  }

  def cleanupFileArray(base: Path) = {
    deleteFile(base.resolve("f1"))
    deleteFile(base.resolve("f2"))
    deleteFile(base.resolve("f3"))
  }

  def deleteFile(path: Path) = Files.delete(path)

  private def createFile(name: String, dir: Path, contents: String) = {
    dir.toFile.mkdirs()
    write(dir.resolve(name).toFile, contents)
  }

  private def write(file: File, contents: String) = {
    val writer = new FileWriter(file)
    writer.write(contents)
    writer.flush()
    writer.close()
    file
  }
}

object SampleWdl {
  object ThreeStep extends SampleWdl {
    override def wdlSource(runtime: String = "") = sourceString()

    private val outputSectionPlaceholder = "OUTPUTSECTIONPLACEHOLDER"

    def sourceString(outputsSection: String = "") = {
      val withPlaceholders =
        """
          |task ps {
          |  command {
          |    ps
          |  }
          |  output {
          |    File procs = stdout()
          |  }
          |}
          |
          |task cgrep {
          |  String pattern
          |  File in_file
          |
          |  command {
          |    grep '${pattern}' ${in_file} | wc -l
          |  }
          |  output {
          |    Int count = read_int(stdout())
          |  }
          |}
          |
          |task wc {
          |  File in_file
          |  command {
          |    cat ${in_file} | wc -l
          |  }
          |  output {
          |    Int count = read_int(stdout())
          |  }
          |}
          |
          |workflow three_step {
          |  call ps
          |  call cgrep {
          |    input: in_file = ps.procs
          |  }
          |  call wc {
          |    input: in_file = ps.procs
          |  }
          |  """ + outputSectionPlaceholder +
          """
          |}
          |
        """
      withPlaceholders.stripMargin.replace(outputSectionPlaceholder, outputsSection)
    }

    val PatternKey = "three_step.cgrep.pattern"
    override lazy val rawInputs = Map(PatternKey -> "...")
  }

  object EmptyInvalid extends SampleWdl {
    override def wdlSource(runtime: String = "") = "{}"
    val rawInputs = Map.empty[String, Any]
  }

  object EmptyWorkflow extends SampleWdl {
    override def wdlSource(runtime: String = "") = "workflow empty_workflow {}"
    val rawInputs = Map.empty[String, Any]
  }

  object EmptyTask extends SampleWdl {
    override def wdlSource(runtime: String = "") = "task empty_task { command { : } }"
    val rawInputs = Map.empty[String, Any]
  }
}
