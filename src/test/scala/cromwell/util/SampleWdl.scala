package cromwell.util

import java.io.{File, FileWriter}
import java.nio.file.{Files, Path}

import wdl4s._
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values._
import cromwell.engine.WorkflowSourceFiles
import spray.json._

import scala.language.postfixOps

trait SampleWdl extends TestFileUtil {
  def wdlSource(runtime: String = ""): WdlSource
  def asWorkflowSources(runtime: String = "", workflowOptions: String = "{}") =
    WorkflowSourceFiles(wdlSource(runtime), wdlJson, workflowOptions)
  val rawInputs: WorkflowRawInputs

  def name = getClass.getSimpleName.stripSuffix("$")

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

  def deleteFile(path: Path) = Files.delete(path)
}

object SampleWdl {

  object HelloWorld extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task hello {
        |  String addressee
        |  command {
        |    echo "Hello ${addressee}!"
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
      """.stripMargin.replaceAll("RUNTIME", runtime)

    val Addressee = "hello.hello.addressee"
    val rawInputs = Map(Addressee -> "world")
    val OutputKey = "hello.hello.salutation"
    val OutputValue = "Hello world!"
  }

  object HelloWorldWithoutWorkflow extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task hello {
        |  String addressee
        |  command {
        |    echo "Hello ${addressee}!"
        |  }
        |  output {
        |    String salutation = read_string(stdout())
        |  }
        |}
      """.stripMargin

    val Addressee = "hello.hello.addressee"
    val rawInputs = Map(Addressee -> "world")
    val OutputKey = "hello.hello.salutation"
    val OutputValue = "Hello world!\n"
  }

  object GoodbyeWorld extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task goodbye {
        |  command {
        |    sh -c "exit 1"
        |  }
        |  output {
        |    String out = read_string(stdout())
        |  }
        |}
        |
        |workflow goodbye {
        |  call goodbye
        |}
      """.stripMargin

    val rawInputs = Map.empty[String, Any]
    val OutputKey = "goodbye.goodbye.out"
  }

  object EmptyWorkflow extends SampleWdl {
    override def wdlSource(runtime: String = "") = "workflow empty_workflow {}"

    val rawInputs = Map.empty[String, Any]
  }

  object EmptyTask extends SampleWdl {
    override def wdlSource(runtime: String = "") = "task empty_task { command { : } }"

    val rawInputs = Map.empty[String, Any]
  }

  object EmptyInvalid extends SampleWdl {
    override def wdlSource(runtime: String = "") = "{}"

    val rawInputs = Map.empty[String, Any]
  }

  object Incr extends SampleWdl {
    override def wdlSource(runtime: String = "") =
    """
      |task incr {
      |  Int val
      |  command {
      |    echo $((${val} + 1))
      |  }
      |  output {
      |    Int out = read_int(stdout())
      |  }
      |}
      |
      |workflow incr {
      |  call incr
      |}
    """.stripMargin

    override val rawInputs: WorkflowRawInputs = Map("incr.incr.val" -> "x1")
  }

  object SubtractionWorkflow {
    val WdlSource =
      """
        |task a {
        |  String message
        |  command { echo '${message}' }
        |  output {
        |    String message = read_string(stdout())
        |    Int constant = 100
        |  }
        |}
        |task b {
        |  command { echo '${message} - ${Int integer}' }
        |}
        |workflow wf {
        |  call a
        |  call b {
        |    input: message = a.message, integer = a.constant - 75
        |  }
        |}
      """.stripMargin
  }

  object CoercionNotDefined extends SampleWdl {
    override def wdlSource(runtime: String = "") = {
      """
        |task summary {
        |  String bfile
        |  command {
        |    ~/plink --bfile ${bfile} --missing --hardy --out foo --allow-no-sex
        |  }
        |  output {
        |    File hwe = "foo.hwe"
        |    File log = "foo.log"
        |    File imiss = "foo.imiss"
        |    File lmiss = "foo.lmiss"
        |  }
        |  meta {
        |    author: "Jackie Goldstein"
        |    email: "jigold@broadinstitute.org"
        |  }
        |}
        |
        |workflow test1 {
        |  call summary {
        |     input: bfile = bfile
        |  }
        |}
      """.stripMargin
    }

    override val rawInputs: WorkflowRawInputs = Map("test1.bfile" -> "data/example1")
  }

  trait ThreeStepTemplate extends SampleWdl {
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
        |  """ + outputSectionPlaceholder + """
        |}
        |
        """
      withPlaceholders.stripMargin.replace(outputSectionPlaceholder, outputsSection)
    }

    val PatternKey ="three_step.cgrep.pattern"
    override lazy val rawInputs = Map(PatternKey -> "...")
  }

  object ThreeStep extends ThreeStepTemplate

  object ThreeStepWithOutputsSection extends ThreeStepTemplate {
    override def wdlSource(runtime: String = "") = sourceString(outputsSection =
      """
        |output {
        | cgrep.count
        | wc.count
        |}
      """.stripMargin)
  }

  object ThreeStepWithInputsInTheOutputsSection extends ThreeStepTemplate {
    override def wdlSource(runtime: String = "") = sourceString(outputsSection =
      """
        |output {
        | cgrep.pattern
        |}
      """.stripMargin)
  }

  object ThreeStepLargeJson extends ThreeStepTemplate {
    override lazy val rawInputs = Map(ThreeStep.PatternKey -> "." * 10000)
  }

  object WorkflowOutputsWithFiles extends SampleWdl {
    // ASCII art from http://www.chris.com/ascii/joan/www.geocities.com/SoHo/7373/flag.html with pipes
    // replaced by exclamation points to keep stripMargin from removing the flagpole.
    override def wdlSource(runtime: String = "") =
      """
        task A {
          command {
            echo "Enfin un peu de francais pour contrer ce raz-de-marée anglais !" > out
            echo "Jacques Chirac fait du jetski sur la Seine en costume traditionnel russe" > out2
          }
          output {
            File out = "out"
            File out2 = "out2"
          }
        }
        task B {
          command {
             echo "Je contre avec un bonnet peruvien et tire une carte chance" > out
             echo "Kamoulox !" > out2
          }
          output {
             Array[File] outs = ["out", "out2"]
          }
        }
        task C {
          command {
            cat > out <<END
            (_)
             !_________________________________________
             !*  *  *  *  * |##########################|
             ! *  *  *  *  *|                          |
             !*  *  *  *  * |##########################|
             ! *  *  *  *  *|                          |
             !*  *  *  *  * |##########################|
             ! *  *  *  *  *|                          |
             !*  *  *  *  * |##########################|
             !~~~~~~~~~~~~~~~                          |
             !#########################################|
             !                                         |
             !#########################################|
             !                                         |
             !###################################JGS###|
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             !
             !
             !
             !
             !
             !
             !
            END
          }
          output {
            File out = "out"
          }
        }
        workflow wfoutputs {
          call A
          call B
          call C
          output {
            A.*
            B.outs
          }
        }
      """.stripMargin
    override lazy val rawInputs = Map.empty[String, String]
  }

  object NestedScatterWdl extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        task A {
        |  command {
        |    echo -n -e "jeff\nchris\nmiguel\nthibault\nkhalid\nscott"
        |  }
        |  output {
        |    Array[String] A_out = read_lines(stdout())
        |  }
        |}
        |
        |task B {
        |  String B_in
        |  command {
        |    python -c "print(len('${B_in}'))"
        |  }
        |  output {
        |    Int B_out = read_int(stdout())
        |  }
        |}
        |
        |task C {
        |  Int C_in
        |  command {
        |    python -c "print(${C_in}*100)"
        |  }
        |  output {
        |    Int C_out = read_int(stdout())
        |  }
        |}
        |
        |task D {
        |  Array[Int] D_in
        |  command {
        |    python -c "print(${sep = '+' D_in})"
        |  }
        |  output {
        |    Int D_out = read_int(stdout())
        |  }
        |}
        |
        |task E {
        |  command {
        |    python -c "import random; print(random.randint(1,100))"
        |  }
        |  output {
        |    Int E_out = read_int(stdout())
        |  }
        |}
        |
        |workflow w {
        |  call A
        |  scatter (item in A.A_out) {
        |    call B {input: B_in = item}
        |    call C {input: C_in = B.B_out}
        |    call E
        |    scatter (itemB in B.B_out) {
        |     call E as G
        |    }
        |    scatter (itemB in B.B_out) {
        |     call E as H
        |    }
        |  }
        |  scatter (item in A.A_out) {
        |    call E as F
        |  }
        |  call D {input: D_in = B.B_out}
        |}
      """.stripMargin
    override lazy val rawInputs = Map("" -> "...")
  }

  object OptionalParamWorkflow extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
    """
      |task hello {
      |  String? person
      |  command {
      |    echo "hello ${person}"
      |  }
      |  output {
      |    String greeting = read_string(stdout())
      |  }
      |}
      |
      |workflow optional {
      |  call hello
      |  call hello as hello2
      |  call hello as hello_person {
      |    input: person = "world"
      |  }
      |}
    """.stripMargin.replaceAll("RUNTIME", runtime)


    override val rawInputs = {
      Map("optional.hello.person" -> "john")
    }
  }

  trait ZeroOrMorePostfixQuantifier extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """
        |task hello {
        |  Array[String] person
        |  command {
        |    echo "hello ${sep = "," person}"
        |  }
        |  output {
        |    String greeting = read_string(stdout())
        |  }
        |}
        |
        |workflow postfix {
        |  call hello
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)
  }

  object ZeroOrMorePostfixQuantifierWorkflowWithArrayInput extends ZeroOrMorePostfixQuantifier {
    override val rawInputs = Map("postfix.hello.person" -> Seq("alice", "bob", "charles"))
  }

  object ZeroOrMorePostfixQuantifierWorkflowWithOneElementArrayInput extends ZeroOrMorePostfixQuantifier {
    override val rawInputs = Map("postfix.hello.person" -> Seq("alice"))
  }

  object ZeroOrMorePostfixQuantifierWorkflowWithZeroElementArrayInput extends ZeroOrMorePostfixQuantifier {
    override val rawInputs = Map("postfix.hello.person" -> Seq())
  }

  trait OneOrMorePostfixQuantifier extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """
        |task hello {
        |  Array[String]+ person
        |  command {
        |    echo "hello ${sep = "," person}"
        |  }
        |  output {
        |    String greeting = read_string(stdout())
        |  }
        |}
        |
        |workflow postfix {
        |  call hello
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)
  }

  object OneOrMorePostfixQuantifierWorkflowWithArrayInput extends OneOrMorePostfixQuantifier {
    override val rawInputs = Map("postfix.hello.person" -> Seq("alice", "bob", "charles"))
  }

  object OneOrMorePostfixQuantifierWorkflowWithScalarInput extends OneOrMorePostfixQuantifier {
    override val rawInputs = Map("postfix.hello.person" -> Seq("alice"))
  }

  trait DefaultParameterValue extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """
        |task hello {
        |  String? person
        |  command {
        |    echo "hello ${default = "default value" person}"
        |  }
        |  output {
        |    String greeting = read_string(stdout())
        |  }
        |}
        |
        |workflow default {
        |  call hello
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)
  }

  object DefaultParameterValueWithValueSpecified extends DefaultParameterValue {
    override val rawInputs = Map("default.hello.person" -> "alice")
  }

  object DefaultParameterValueWithNOValueSpecified extends DefaultParameterValue {
    override val rawInputs = Map.empty[String, String]
  }
  object CannedThreeStep extends SampleWdl {
    val CannedProcessOutput =
      """
        |USER              PID  %CPU %MEM      VSZ    RSS   TT  STAT STARTED      TIME COMMAND
        |joeblaux         10344   4.5 20.6  7011548 3454616  ??  S    Mon06AM 275:26.10 /Applications/IntelliJ IDEA 14.app/Contents/MacOS/idea
        |joeblaux           883   2.2  0.5  2716336  85768   ??  S    Sun08AM   0:52.64 /Applications/Utilities/Terminal.app/Contents/MacOS/Terminal
        |_coreaudiod        577   1.9  0.1  2522428   9572   ??  Ss   Sun08AM  55:39.69 /usr/sbin/coreaudiod
        |_windowserver      130   1.6  1.9  4535136 319588   ??  Ss   Sun08AM 148:39.39 /System/Library/Frameworks/ApplicationServices.framework/Frameworks/CoreGraphics.framework/Resources/WindowServer -daemon
        |joeblaux          6440   1.4  2.2  4496164 362136   ??  S    Sun09PM  74:29.40 /Applications/Google Chrome.app/Contents/MacOS/Google Chrome
      """.stripMargin.trim

    override def wdlSource(runtime: String): WdlSource =
      """
        |task ps {
        |  File dummy_ps_file
        |  command {
        |    cat ${dummy_ps_file}
        |  }
        |  output {
        |    File procs = stdout()
        |  }
        |  RUNTIME
        |}
        |
        |task cgrep {
        |  String pattern
        |  File in_file
        |  command {
        |    grep '${pattern}' ${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |  RUNTIME
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
        |  RUNTIME
        |}
        |
        |workflow three_step {
        |  call ps
        |  call ps as ps2
        |  call ps as ps3
        |  call cgrep {
        |    input: in_file = ps.procs
        |  }
        |  call wc {
        |    input: in_file = ps.procs
        |  }
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    private def createCannedPsFile: File = {
      val file = File.createTempFile("canned_ps", ".out")
      val writer = new FileWriter(file)
      writer.write(CannedProcessOutput)
      writer.flush()
      writer.close()
      file
    }

    override val rawInputs = {
      Map(
        "three_step.cgrep.pattern" -> "joeblaux",
        "three_step.ps.dummy_ps_file" -> createCannedPsFile.getAbsolutePath,
        "three_step.ps2.dummy_ps_file" -> createCannedPsFile.getAbsolutePath,
        "three_step.ps3.dummy_ps_file" -> createCannedPsFile.getAbsolutePath
      )
    }
  }

  object CannedFilePassing extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task ps {
        |  command {
        |    (echo "x"; echo "y"; echo "z") > myfile.txt
        |  }
        |  output {
        |    File procs = "myfile.txt"
        |  }
        |  runtime {
        |    docker: "ubuntu:latest"
        |  }
        |}
        |
        |task cgrep {
        |  String pattern
        |  File in_file
        |  command {
        |    grep '${pattern}' ${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |  runtime {
        |    docker: "ubuntu:latest"
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
        |  runtime {
        |    docker: "ubuntu:latest"
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
        |}
        |
      """.stripMargin

    override val rawInputs = Map(ThreeStep.PatternKey -> "x")
  }

  object CurrentDirectory extends SampleWdl {
    override def wdlSource(runtime: String): String =
      """
        |task whereami {
        |  command {
        |    pwd
        |  }
        |  output {
        |    String pwd = read_string(stdout())
        |  }
        |}
        |
        |workflow whereami {
        |  call whereami
        |}
      """.stripMargin

    override val rawInputs: Map[String, Any] = Map.empty
  }

  object ReadLinesFunctionWdl extends SampleWdl {
    val CannedOutput =
      """java
        |scala
        |c
        |c++
        |python
        |bash
      """.stripMargin.trim

    override def wdlSource(runtime: String): WdlSource =
      """
        |task cat_to_stdout {
        |  File file
        |  command {
        |    cat ${file}
        |  }
        |  output {
        |    Array[String] lines = read_lines(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |task cat_to_file {
        |  File file
        |  command {
        |    cat ${file} > out
        |  }
        |  output {
        |    Array[String] lines = read_lines("out")
        |  }
        |  RUNTIME
        |}
        |
        |workflow read_lines {
        |  call cat_to_stdout
        |  call cat_to_file
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    override val rawInputs = {
      Map(
        "read_lines.cat_to_stdout.file" -> createCannedFile("canned", CannedOutput).getAbsolutePath,
        "read_lines.cat_to_file.file" -> createCannedFile("canned1", CannedOutput).getAbsolutePath
      )
    }
  }

  object DeclarationsWorkflow extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """
        |task cat {
        |  File file
        |  String? flags
        |  command {
        |    cat ${flags} ${file}
        |  }
        |  output {
        |    File procs = stdout()
        |  }
        |}
        |
        |task cgrep {
        |  String str_decl
        |  String pattern
        |  File in_file
        |  command {
        |    grep '${pattern}' ${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |    String str = str_decl
        |  }
        |}
        |
        |workflow two_step {
        |  String flags_suffix
        |  String flags = "-" + flags_suffix
        |  String static_string = "foobarbaz"
        |  call cat {
        |    input: flags = flags
        |  }
        |  call cgrep {
        |    input: in_file = cat.procs
        |  }
        |}
      """.stripMargin

    private val fileContents =
      s"""first line
        |second line
        |third line
       """.stripMargin

    override val rawInputs: WorkflowRawInputs = Map(
      "two_step.cgrep.pattern" -> "first",
      "two_step.cgrep.str_decl" -> "foobar",
      "two_step.cat.file" -> createCannedFile("canned", fileContents).getAbsolutePath,
      "two_step.flags_suffix" -> "s"
    )
  }

  object StringInterpolation extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task echo {
        |  String greeting
        |  String out
        |  Int one = 1
        |
        |  command {
        |    echo "${greeting}" > ${out}.txt
        |    echo "${ one + 1 }" > ${one+1}.txt
        |  }
        |  output {
        |    File outfile = "${ out }.txt"
        |    Int two = read_int("${ one + 1 }.txt")
        |  }
        |}
        |
        |workflow echo_wf {
        |  call echo
        |}
        |
      """.stripMargin

    override val rawInputs = Map(
      "echo_wf.echo.greeting" -> "world",
      "echo_wf.echo.out" -> "foobar"
    )
  }

  object ArrayIO extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task concat_files {
        |  String? flags
        |  Array[File]+ files
        |  command {
        |    cat ${default = "-s" flags} ${sep = " " files}
        |  }
        |  output {
        |    File concatenated = stdout()
        |  }
        |}
        |
        |task find {
        |  String pattern
        |  File root
        |  command {
        |    find ${root} ${"-name " + pattern}
        |  }
        |  output {
        |    Array[String] results = read_lines(stdout())
        |  }
        |}
        |
        |task count_lines {
        |  Array[File]+ files
        |  command {
        |    cat ${sep = ' ' files} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |}
        |
        |task serialize {
        |  Array[String] strs
        |  command {
        |    cat ${write_lines(strs)}
        |  }
        |  output {
        |    String contents = read_string(stdout())
        |  }
        |}
        |
        |workflow wf {
        |  Array[File] files
        |  Array[String] strings = ["str1", "str2", "str3"]
        |  call serialize {
        |    input: strs = strings
        |  }
        |  call concat_files as concat {
        |    input: files = files
        |  }
        |  call count_lines {
        |    input: files = [concat.concatenated]
        |  }
        |  call find
        |  call count_lines as count_lines_array {
        |    input: files = find.results
        |  }
        |}
      """.stripMargin

    val tempDir = Files.createTempDirectory("ArrayIO")
    val firstFile = createCannedFile(prefix = "first", contents = "foo\n", dir = Some(tempDir))
    val secondFile = createCannedFile(prefix = "second", contents = "bar\nbaz\n", dir = Some(tempDir))

    override val rawInputs = Map(
      "wf.find.root" -> tempDir.toAbsolutePath.toString,
      "wf.find.pattern" -> "*.out", // createCannedFile makes files that have .out extension
      "wf.files" -> Seq(firstFile.getAbsolutePath, secondFile.getAbsolutePath)
    )
  }

  case class ArrayLiteral(catRootDir: Path) extends SampleWdl {
    createFileArray(catRootDir)
    def cleanup() = cleanupFileArray(catRootDir)

    override def wdlSource(runtime: String = "") =
      """
        |task cat {
        |  Array[File]+ files
        |  command {
        |    cat -s ${sep = ' ' files}
        |  }
        |  output {
        |    Array[String] lines = read_lines(stdout())
        |  }
        |}
        |
        |workflow wf {
        |  Array[File] arr = ["f1", "f2", "f3"]
        |  call cat {input: files = arr}
        |}
      """.stripMargin

    override val rawInputs = Map.empty[String, String]
  }

  case class MapLiteral(catRootDir: Path) extends SampleWdl {
    createFileArray(catRootDir)
    def cleanup() = cleanupFileArray(catRootDir)

    override def wdlSource(runtime: String = "") =
      """
        |task write_map {
        |  Map[File, String] file_to_name
        |  command {
        |    cat ${write_map(file_to_name)}
        |  }
        |  output {
        |    String contents = read_string(stdout())
        |  }
        |}
        |
        |task read_map {
        |  command <<<
        |    python <<CODE
        |    map = {'x': 500, 'y': 600, 'z': 700}
        |    print("\n".join(["{}\t{}".format(k,v) for k,v in map.items()]))
        |    CODE
        |  >>>
        |  output {
        |    Map[String, Int] out_map = read_map(stdout())
        |  }
        |}
        |
        |workflow wf {
        |  Map[File, String] map = {"f1": "alice", "f2": "bob", "f3": "chuck"}
        |  call write_map {input: file_to_name = map}
        |  call read_map
        |}
      """.stripMargin

    override val rawInputs = Map.empty[String, String]
  }

  object MultiLineCommandWorkflowWdl extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task blah {
        |  command <<<
        |    python <<CODE
        |    def a():
        |      return "a"
        |    def b():
        |      return "b"
        |    print('{}{}'.format(a(),b()))
        |    CODE
        |  >>>
        |
        |  output {
        |    String ab = read_string(stdout())
        |  }
        |}
        |
        |workflow wf {
        |  call blah
        |}
      """.stripMargin

    override val rawInputs = Map.empty[String, String]
  }

  object InputIsolationWdl extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task localize {
        |    Array[File] array
        |    command {
        |        ls -1 "$(dirname ${array[0]})" | wc -l | tr -d '[[:space:]]'
        |    }
        |    output {
        |        String ls = read_string(stdout())
        |    }
        |    RUNTIME
        |}
        |
        |task echo_int {
        |  Int int
        |  command {echo ${int} > out }
        |  output {File out = "out"}
        |   RUNTIME
        |}
        |
        |workflow wf {
        |    Array[File] files
        |    Array[Int] ints = [1,2]
        |
        |   scatter(i in ints) {
        |    call echo_int {
        |      input: int = i
        |    }
        |  }
        |
        |  call localize as fromDifferentDirectories { input: array = echo_int.out }
        |  call localize as fromSameDirectory { input: array = files }
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    val tempDir = Files.createTempDirectory("InputFiles")
    val file1 = createCannedFile(prefix = "file1", contents = "", dir = Some(tempDir))
    val file2 = createCannedFile(prefix = "file1", contents = "", dir = Some(tempDir))

    override val rawInputs = Map(
      "wf.files" -> Seq(file1.getAbsolutePath, file2.getAbsolutePath)
    )
  }

  object TripleSleep extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task wait {
        |  command <<<
        |    sleep 100
        |    echo "waited 100 seconds"
        |  >>>
        |  output {
        |    String intermediate = read_string(stdout())
        |  }
        |}
        |task msg {
        |  command {
        |    echo ${sommat}
        |  }
        |  output {
        |    String result = read_string(stdout())
        |  }
        |}
        |
        |workflow hello {
        |  call wait as wait1
        |  call wait as wait2
        |  call wait as wait3
        |
        |  call msg as msg1 {input: sommat = wait1.intermediate}
        |  call msg as msg2 {input: sommat = wait2.intermediate}
        |  call msg as msg3 {input: sommat = wait3.intermediate}
        |}
      """.stripMargin

    override val rawInputs = Map.empty[String, String]
    val OutputKey1 = "hello.msg1.result"
    val OutputValue = "waited 100 seconds"
  }

  object BadTaskOutputWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task bad {
        |  command {
        |    echo "hello" > a
        |  }
        |  output {
        |    # Oops! we made a spelling mistake in our WDL!
        |    File a = "b"
        |  }
        |}
        |
        |workflow badExample {
        |  call bad
        |}
      """.stripMargin

    override val rawInputs =  Map.empty[String, String]
  }

  object ArrayAndMapIndexingWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task echo_str {
        |  String s
        |  command { echo ${s} }
        |  output { String o = read_string(stdout()) }
        |  RUNTIME
        |}
        |
        |task echo_int {
        |  Int i
        |  command { echo ${i} }
        |  output { Int o = read_int(stdout()) }
        |  RUNTIME
        |}
        |
        |workflow test {
        |  Map[String, Int] m = {"a": 0, "b": 1, "c": 2}
        |  Array[String] a = ["foo", "bar", "baz"]
        |  call echo_str {
        |    input: s = a[1]
        |  }
        |  call echo_int {
        |    input: i = m["c"]*100
        |  }
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    override val rawInputs =  Map.empty[String, String]
  }

  class ScatterWdl extends SampleWdl {
    val tasks = """task A {
      |  command {
      |    echo -n -e "jeff\nchris\nmiguel\nthibault\nkhalid\nscott"
      |  }
      |  RUNTIME
      |  output {
      |    Array[String] A_out = read_lines(stdout())
      |  }
      |}
      |
      |task B {
      |  String B_in
      |  command {
      |    python -c "print(len('${B_in}'))"
      |  }
      |  RUNTIME
      |  output {
      |    Int B_out = read_int(stdout())
      |  }
      |}
      |
      |task C {
      |  Int C_in
      |  command {
      |    python -c "print(${C_in}*100)"
      |  }
      |  RUNTIME
      |  output {
      |    Int C_out = read_int(stdout())
      |  }
      |}
      |
      |task D {
      |  Array[Int] D_in
      |  command {
      |    python -c "print(${sep = '+' D_in})"
      |  }
      |  RUNTIME
      |  output {
      |    Int D_out = read_int(stdout())
      |  }
      |}
      |
      |task E {
      |  command {
      |    python -c "print(9)"
      |  }
      |  RUNTIME
      |  output {
      |    Int E_out = read_int(stdout())
      |  }
      |}
    """.stripMargin

    override def wdlSource(runtime: String = "") =
      s"""$tasks
        |
        |workflow w {
        |  call A
        |  scatter (item in A.A_out) {
        |    call B {input: B_in = item}
        |    call C {input: C_in = B.B_out}
        |    call E
        |  }
        |  call D {input: D_in = B.B_out}
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    override lazy val rawInputs = Map.empty[String, String]
  }

  object SiblingsScatterWdl extends ScatterWdl {
    override def wdlSource(runtime: String = "") =
      s"""$tasks
        |
        |workflow w {
        |  call A
        |  scatter (item in A.A_out) {
        |    call B {input: B_in = item}
        |    call C {input: C_in = B.B_out}
        |    call E
        |  }
        |  scatter (item in A.A_out) {
        |    call B as F {input: B_in = item}
        |  }
        |  call D {input: D_in = B.B_out}
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    override lazy val rawInputs = Map.empty[String, String]
  }

  object SimpleScatterWdl extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task echo_int {
        |  Int int
        |  command {echo ${int}}
        |  output {Int out = read_int(stdout())}
        |}
        |
        |workflow scatter0 {
        |  Array[Int] ints = [1,2,3,4,5]
        |  call echo_int as outside_scatter {input: int = 8000}
        |  scatter(i in ints) {
        |    call echo_int as inside_scatter {
        |      input: int = i
        |    }
        |  }
        |}
      """.stripMargin

    override lazy val rawInputs = Map.empty[String, String]
  }

  object SimpleScatterWdlWithOutputs extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task echo_int {
        |  Int int
        |  command {echo ${int}}
        |  output {Int out = read_int(stdout())}
        |}
        |
        |workflow scatter0 {
        |  Array[Int] ints = [1,2,3,4,5]
        |  call echo_int as outside_scatter {input: int = 8000}
        |  scatter(i in ints) {
        |    call echo_int as inside_scatter {
        |      input: int = i
        |    }
        |  }
        |  output {
        |    inside_scatter.*
        |  }
        |}
      """.stripMargin

    override lazy val rawInputs = Map.empty[String, String]
  }

  object PrepareScatterGatherWdl extends SampleWdl {
    override def wdlSource(runtime: String = "") = {
      """
        |#
        |# Goal here is to split up the input file into files of 1 line each (in the prepare) then in parallel call wc -w on each newly created file and count the words into another file then in the gather, sum the results of each parallel call to come up with
        |# the word-count for the fil
        |#
        |# splits each line into a file with the name temp_?? (shuffle)
        |task do_prepare {
        |    File input_file
        |    command {
        |        split -l 1 ${input_file} temp_ && ls -1 temp_?? > files.list
        |    }
        |    output {
        |        Array[File] split_files = read_lines("files.list")
        |    }
        |    RUNTIME
        |}
        |# count the number of words in the input file, writing the count to an output file overkill in this case, but simulates a real scatter-gather that would just return an Int (map)
        |task do_scatter {
        |    File input_file
        |    command {
        |        wc -w ${input_file} > output.txt
        |    }
        |    output {
        |        File count_file = "output.txt"
        |    }
        |    RUNTIME
        |}
        |# aggregate the results back together (reduce)
        |task do_gather {
        |    Array[File] input_files
        |    command <<<
        |        cat ${sep = ' ' input_files} | awk '{s+=$1} END {print s}'
        |    >>>
        |    output {
        |        Int sum = read_int(stdout())
        |    }
        |    RUNTIME
        |}
        |workflow sc_test {
        |    call do_prepare
        |    scatter(f in do_prepare.split_files) {
        |        call do_scatter {
        |            input: input_file = f
        |        }
        |    }
        |    call do_gather {
        |        input: input_files = do_scatter.count_file
        |    }
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)
    }

    val contents =
        """|the
           |total number
           |of words in this
           |text file is 11
           |""".stripMargin

    override lazy val rawInputs = Map("sc_test.do_prepare.input_file" -> createCannedFile("scatter",contents).getAbsolutePath)
  }

  object FileClobber extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task read_line {
        |  File in
        |  command { cat ${in} }
        |  output { String out = read_string(stdout()) }
        |}
        |
        |workflow two {
        |  call read_line as x
        |  call read_line as y
        |}
      """.stripMargin

    val tempDir1 = Files.createTempDirectory("FileClobber1")
    val tempDir2 = Files.createTempDirectory("FileClobber2")
    val firstFile = createFile(name = "file.txt", contents = "first file.txt", dir = tempDir1)
    val secondFile = createFile(name = "file.txt", contents = "second file.txt", dir = tempDir2)

    override val rawInputs = Map(
      "two.x.in" -> firstFile.getAbsolutePath,
      "two.y.in" -> secondFile.getAbsolutePath
    )
  }

  object ContinueOnReturnCode extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        task A {
        |  command {
        |    python -c "print(321);exit(123)"
        |  }
        |  output {
        |    Int A_out = read_int(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |task B {
        |  Int B_in
        |  command {
        |    echo ${B_in}
        |  }
        |}
        |
        |
        |workflow w {
        |  call A
        |  call B {input: B_in = A.A_out}
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)
    override lazy val rawInputs = Map("" -> "...")
  }

  object FilePassingWorkflow extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task a {
        |  File in
        |  String out_name = "out"
        |
        |  command {
        |    cat ${in} > ${out_name}
        |  }
        |  RUNTIME
        |  output {
        |    File out = "out"
        |    File out_interpolation = "${out_name}"
        |    String contents = read_string("${out_name}")
        |  }
        |}
        |
        |workflow file_passing {
        |  File f
        |
        |  call a {input: in = f}
        |  call a as b {input: in = a.out}
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    private val fileContents = s"foo bar baz"

    override val rawInputs: WorkflowRawInputs = Map(
      "file_passing.f" -> createCannedFile("canned", fileContents).getAbsolutePath
    )
  }

  /**
    * @param salt - an arbitrary value that will be added as
    *               a BASH comment on the command, this is so
    *               tests can have control over call caching
    *               for this workflow.  i.e. so one test can't
    *               call cache to another test if the seeds are
    *               different
    */
  case class CallCachingWorkflow(salt: String) extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task a {
        |  File in
        |  String out_name = "out"
        |  String salt
        |
        |  command {
        |    # ${salt}
        |    cat ${in} > ${out_name}
        |  }
        |  RUNTIME
        |  output {
        |    File out = "out"
        |    File out_interpolation = "${out_name}"
        |    String contents = read_string("${out_name}")
        |  }
        |}
        |
        |workflow file_passing {
        |  File f
        |
        |  call a {input: in = f}
        |  call a as b {input: in = a.out}
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    private val fileContents = s"foo bar baz"

    override val rawInputs: WorkflowRawInputs = Map(
      "file_passing.f" -> createCannedFile("canned", fileContents).getAbsolutePath,
      "file_passing.a.salt" -> salt,
      "file_passing.b.salt" -> salt
    )
  }

  object WdlFunctionsAtWorkflowLevel extends SampleWdl {
    val CannedArray =
      """one
        |two
        |three
        |four
        |five
      """.stripMargin.trim

    val CannedMap =
      s"""k1\tv1
        |k2\tv2
        |k3\tv3
      """.stripMargin.trim

    override def wdlSource(runtime: String): WdlSource =
      """
        |task a {
        |  Array[String] array
        |  Map[String, String] map
        |
        |  command {
        |    echo ${sep = ' ' array} > concat
        |  }
        |  output {
        |    String x = read_string("concat")
        |    Map[String, String] y = map
        |  }
        |}
        |
        |workflow w {
        |  File array_file
        |  File map_file
        |  Array[String] in_array = read_lines(array_file)
        |  Map[String, String] in_map = read_map(map_file)
        |  call a {input:
        |    array = in_array,
        |    map = in_map
        |  }
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    override val rawInputs = {
      Map(
        "w.array_file" -> createCannedFile("array.txt", CannedArray).getAbsolutePath,
        "w.map_file" -> createCannedFile("map.txt", CannedMap).getAbsolutePath
      )
    }
  }

  object ArrayOfArrays extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task subtask {
        |  Array[File] a
        |  command {
        |    cat ${sep = " " a}
        |  }
        |  output {
        |    String concatenated = read_string(stdout())
        |  }
        |}
        |
        |workflow wf {
        |  Array[Array[File]] nested_file
        |
        |  scatter(n in nested_file) {
        |    call subtask {
        |      input: a = n
        |    }
        |  }
        |}
      """.stripMargin

    val tempDir = Files.createTempDirectory("ArrayOfArray")
    val firstFile = createCannedFile(prefix = "first", contents = "foo\n", dir = Some(tempDir))
    val secondFile = createCannedFile(prefix = "second", contents = "bar\nbaz\n", dir = Some(tempDir))
    val thirdFile = createCannedFile(prefix = "third", contents = "third\n", dir = Some(tempDir))
    val fourthFile = createCannedFile(prefix = "fourth", contents = "fourth\n", dir = Some(tempDir))

    override val rawInputs = Map(
      "wf.nested_file" ->
        WdlArray(WdlArrayType(WdlArrayType(WdlStringType)),
        Seq(
          WdlArray(WdlArrayType(WdlStringType), Seq(firstFile.getAbsolutePath, secondFile.getAbsolutePath).map(WdlString)),
          WdlArray(WdlArrayType(WdlStringType), Seq(thirdFile.getAbsolutePath, fourthFile.getAbsolutePath).map(WdlString))
        )
      )
    )
  }

  object GlobtasticWorkflow extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        task A {
        |  command {
        |    mkdir out
        |    echo a > out/a.txt
        |    echo b > out/b.txt
        |    echo c > out/c.txt
        |  }
        |  output {
        |    Array[File] A_out = glob("out/*.txt")
        |  }
        |  RUNTIME
        |}
        |
        |task B {
        |  Array[File] B_in
        |  command {
        |    cat ${sep=' ' B_in}
        |  }
        |  output {
        |    String B_out = read_string(stdout())
        |  }
        |}
        |
        |
        |workflow w {
        |  call A
        |  call B {input: B_in=A.A_out}
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)
    override lazy val rawInputs = Map("" -> "...")
  }

  object ReadTsvWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task output_table {
        |  command {
        |    echo -e "col0\tcol1\tcol2"
        |    echo -e "a\tb\tc"
        |    echo -e "x\ty\tz"
        |  }
        |  output {
        |     Array[Array[String]] table = read_tsv(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |task output_file_table {
        |  command {
        |    echo "first" > first
        |    echo "second" > second
        |    echo "third" > third
        |    echo "fourth" > fourth
        |    echo -e "first\tsecond"
        |    echo -e "third\tfourth"
        |  }
        |  output {
        |     Array[Array[File]] table = read_tsv(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |task output_matrix {
        |  command {
        |    echo -e "0\t1\t2"
        |    echo -e "3\t4\t5"
        |    echo -e "6\t7\t8"
        |  }
        |  output {
        |     Array[Array[Int]] matrix = read_tsv(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |workflow test {
        |  call output_table
        |  call output_file_table
        |  call output_matrix
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    override val rawInputs = Map.empty[String, String]
  }

  object CallCachingHashingWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task t {
        |  Int a
        |  Float b
        |  String c
        |  File d
        |
        |  command {
        |    echo "${a}" > a
        |    echo "${b}" > b
        |    echo "${c}" > c
        |    cat ${d} > d
        |  }
        |  output {
        |    Int w = read_int("a") + 2
        |    Float x = read_float("b")
        |    String y = read_string("c")
        |    File z = "d"
        |  }
        |  RUNTIME
        |}
        |
        |workflow w {
        |  call t
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    val tempDir = Files.createTempDirectory("CallCachingHashingWdl")
    val cannedFile = createCannedFile(prefix = "canned", contents = "file contents", dir = Some(tempDir))
    override val rawInputs = Map(
      "w.t.a" -> WdlInteger(1),
      "w.t.b" -> WdlFloat(1.1),
      "w.t.c" -> WdlString("foobar"),
      "w.t.d" -> WdlFile(cannedFile.getAbsolutePath)
    )
  }

  object ExpressionsInInputs extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task echo {
        |  String inString
        |  command {
        |    echo ${inString}
        |  }
        |
        |  output {
        |    String outString = read_string(stdout())
        |  }
        |}
        |
        |workflow wf {
        |  String a1
        |  String a2
        |  call echo {
        |   input: inString = a1 + " " + a2
        |  }
        |  call echo as echo2 {
        |    input: inString = a1 + " " + echo.outString + " " + a2
        |  }
        |}
      """.stripMargin
    override val rawInputs = Map(
      "wf.a1" -> WdlString("hello"),
      "wf.a2" -> WdlString("world")
    )
  }

  object SingleToArrayCoercion extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task singleFile {
        |  command {
        |    echo hello
        |  }
        |  output {
        |    File out = stdout()
        |  }
        |}
        |
        |task listFiles {
        |  Array[File] manyIn
        |  command {
        |    cat ${sep=" " manyIn}
        |  }
        |  output {
        |    String result = read_string(stdout())
        |  }
        |}
        |
        |workflow oneToMany {
        |  call singleFile
        |  call listFiles { input: manyIn = singleFile.out }
        |}
      """.stripMargin
    override val rawInputs = Map.empty[String, String]
  }

  /**
    * Inputs referencing other inputs and outputs referencing other outputs.
    */
  object ReferencingPreviousInputsAndOutputs extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task golden_pie {
        |  Float pi = 3.1415926
        |  Float tau = pi + pi
        |
        |  command {
        |    echo 1.6180339887
        |    echo ${tau} 1>&2
        |  }
        |
        |  output {
        |    Float Au = read_float(stdout())
        |    Float doubleAu = Au + Au
        |    Float tauValue = read_float(stderr())
        |    Float goldenPie = Au * pi
        |  }
        |}
        |
        |workflow wf {
        |  call golden_pie
        |}
      """.stripMargin
    override val rawInputs = Map.empty[FullyQualifiedName, Any]
  }
}
