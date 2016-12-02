package cromwell.util

import java.nio.file.{Files, Path}
import java.util.UUID

import better.files._
import cromwell.core.{WorkflowSourceFilesWithoutImports}
import spray.json._
import wdl4s._
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values._

import scala.language.postfixOps

trait SampleWdl extends TestFileUtil {
  def wdlSource(runtime: String = ""): WdlSource
  def asWorkflowSources(runtime: String = "", workflowOptions: String = "{}") =
    WorkflowSourceFilesWithoutImports(wdlSource(runtime), wdlJson, workflowOptions)
  val rawInputs: WorkflowRawInputs

  def name = getClass.getSimpleName.stripSuffix("$")

  def createFileArray(base: Path): Unit = {
    createFile("f1", base, "line1\nline2\n")
    createFile("f2", base, "line3\nline4\n")
    createFile("f3", base, "line5\n")
    ()
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
      s"""
        |task hello {
        |  String addressee
        |  command {
        |    echo "Hello $${addressee}!"
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
      """.stripMargin.replaceAll("RUNTIME", runtime)

    val Addressee = "wf_hello.hello.addressee"
    val rawInputs = Map(Addressee -> "world")
    val OutputKey = "wf_hello.hello.salutation"
    val OutputValue = "Hello world!"
  }

  object HelloWorldWithoutWorkflow extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      s"""
        |task hello {
        |  String addressee
        |  command {
        |    echo "Hello $${addressee}!"
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
        |workflow wf_goodbye {
        |  call goodbye
        |}
      """.stripMargin

    val rawInputs = Map.empty[String, Any]
    val OutputKey = "goodbye.goodbye.out"
  }

  object EmptyString extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      s"""
        |task hello {
        |  command {
        |    echo "Hello!"
        |  }
        |  output {
        |    String empty = ""
        |  }
        |  RUNTIME
        |}
        |
        |task goodbye {
        |  String emptyInputString
        |  command {
        |    echo "$${emptyInputString}"
        |  }
        |  output {
        |    String empty = read_string(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |workflow wf_hello {
        |  call hello
        |  call goodbye {input: emptyInputString=hello.empty }
        |  output {
        |   hello.empty
        |   goodbye.empty
        |  }
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    val rawInputs = Map.empty[String, Any]
    val outputMap = Map(
      "hello.hello.empty" -> WdlString(""),
      "hello.goodbye.empty" -> WdlString("")
    )
  }


  object EmptyWorkflow extends SampleWdl {
    override def wdlSource(runtime: String = "") = "workflow empty_workflow {}"

    val rawInputs = Map.empty[String, Any]
  }

  object CoercionNotDefined extends SampleWdl {
    override def wdlSource(runtime: String = "") = {
      s"""
        |task summary {
        |  String bfile
        |  command {
        |    ~/plink --bfile $${bfile} --missing --hardy --out foo --allow-no-sex
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
    override def wdlSource(runtime: String = "") = sourceString().replaceAll("RUNTIME", runtime)
    private val outputSectionPlaceholder = "OUTPUTSECTIONPLACEHOLDER"
    def sourceString(outputsSection: String = "") = {
      val withPlaceholders =
        s"""
        |task ps {
        |  command {
        |    ps
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
        |
        |  command {
        |    grep '$${pattern}' $${in_file} | wc -l
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
        |    cat $${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |  RUNTIME
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
      """.stripMargin).replaceAll("RUNTIME", runtime)
  }

  object ThreeStepWithInputsInTheOutputsSection extends ThreeStepTemplate {
    override def wdlSource(runtime: String = "") = sourceString(outputsSection =
      """
        |output {
        | cgrep.pattern
        |}
      """.stripMargin).replaceAll("RUNTIME", runtime)
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

  object WorkflowScatterOutputsWithFileArrays extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task A {
        |  command {
        |    # Make sure that even for slow tasks, the output-copy waits until the tasks are complete before starting to copy.
        |    sleep 0.1
        |    echo "The creatures outside looked from pig to man, and from man to pig, and from pig to man again: but already it was impossible to say which was which." > B1
        |    echo "But it was all right, everything was all right, the struggle was finished. He had won the victory over himself. He loved Big Brother." > B2
        |  }
        |  output {
        |    Array[File] outs = [ "B1", "B2" ]
        |  }
        |}
        |
        |workflow wfoutputs {
        |  Array[Int] xs = [1,2,3,4,5,6,7,8,9,10]
        |
        |  scatter ( x in xs ) {
        |    call A
        |  }
        |}
      """.stripMargin
    override lazy val rawInputs = Map.empty[String, String]
  }


  object DeclarationsWorkflow extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      s"""
        |task cat {
        |  File file
        |  String? flags
        |  command {
        |    cat $${flags} $${file}
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
        |    grep '$${pattern}' $${in_file} | wc -l
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
      "two_step.cat.file" -> createCannedFile("canned", fileContents),
      "two_step.flags_suffix" -> "s"
    )
  }

  trait ZeroOrMorePostfixQuantifier extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      s"""
        |task hello {
        |  Array[String] person
        |  command {
        |    echo "hello $${sep = "," person}"
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
      s"""
        |task hello {
        |  Array[String]+ person
        |  command {
        |    echo "hello $${sep = "," person}"
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
        |  RUNTIME
        |}
        |
        |workflow wf_whereami {
        |  call whereami
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)

    override val rawInputs: Map[String, Any] = Map.empty
  }

  object ArrayIO extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      s"""task concat_files {
        |  String? flags
        |  Array[File]+ files
        |  command {
        |    cat $${default = "-s" flags} $${sep = " " files}
        |  }
        |  output {
        |    File concatenated = stdout()
        |  }
        |  RUNTIME
        |}
        |
        |task find {
        |  String pattern
        |  File root
        |  command {
        |    find $${root} $${"-name " + pattern}
        |  }
        |  output {
        |    Array[String] results = read_lines(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |task count_lines {
        |  Array[File]+ files
        |  command {
        |    cat $${sep = ' ' files} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |task serialize {
        |  Array[String] strs
        |  command {
        |    cat $${write_lines(strs)}
        |  }
        |  output {
        |    String contents = read_string(stdout())
        |  }
        |  RUNTIME
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
      """.stripMargin.replaceAll("RUNTIME", runtime)

    val tempDir = File.newTemporaryDirectory("ArrayIO").path
    val firstFile = createCannedFile(prefix = "first", contents = "foo\n", dir = Some(tempDir))
    val secondFile = createCannedFile(prefix = "second", contents = "bar\nbaz\n", dir = Some(tempDir))

    override val rawInputs = Map(
      "wf.find.root" -> tempDir.toString,
      "wf.find.pattern" -> "*.out", // createCannedFile makes files that have .out extension
      "wf.files" -> Seq(firstFile.toString, secondFile.toString)
    )
  }

  case class ArrayLiteral(catRootDir: Path) extends SampleWdl {
    createFileArray(catRootDir)
    def cleanup() = cleanupFileArray(catRootDir)

    override def wdlSource(runtime: String = "") =
      s"""
        |task cat {
        |  Array[File]+ files
        |  command {
        |    cat -s $${sep = ' ' files}
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
      s"""
        |task write_map {
        |  Map[File, String] file_to_name
        |  command {
        |    cat $${write_map(file_to_name)}
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
        |    print("\\n".join(["{}\\t{}".format(k,v) for k,v in map.items()]))
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

  class ScatterWdl extends SampleWdl {
    val tasks = s"""task A {
      |  command {
      |    echo -n -e "jeff\nchris\nmiguel\nthibault\nkhalid\nruchi"
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
      |    python -c "print(len('$${B_in}'))"
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
      |    python -c "print($${C_in}*100)"
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
      |    python -c "print($${sep = '+' D_in})"
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
      s"""task echo_int {
        |  Int int
        |  command {echo $${int}}
        |  output {Int out = read_int(stdout())}
        |  RUNTIME_PLACEHOLDER
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
      """.stripMargin.replace("RUNTIME_PLACEHOLDER", runtime)

    override lazy val rawInputs = Map.empty[String, String]
  }

  object SimpleScatterWdlWithOutputs extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      s"""task echo_int {
        |  Int int
        |  command {echo $${int}}
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

  case class PrepareScatterGatherWdl(salt: String = UUID.randomUUID().toString) extends SampleWdl {
    override def wdlSource(runtime: String = "") = {
      s"""
        |#
        |# Goal here is to split up the input file into files of 1 line each (in the prepare) then in parallel call wc -w on each newly created file and count the words into another file then in the gather, sum the results of each parallel call to come up with
        |# the word-count for the fil
        |#
        |# splits each line into a file with the name temp_?? (shuffle)
        |task do_prepare {
        |    File input_file
        |    command {
        |        split -l 1 $${input_file} temp_ && ls -1 temp_?? > files.list
        |    }
        |    output {
        |        Array[File] split_files = read_lines("files.list")
        |    }
        |    RUNTIME
        |}
        |# count the number of words in the input file, writing the count to an output file overkill in this case, but simulates a real scatter-gather that would just return an Int (map)
        |task do_scatter {
        |    String salt
        |    File input_file
        |    command {
        |        # $${salt}
        |        wc -w $${input_file} > output.txt
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
        |        cat $${sep = ' ' input_files} | awk '{s+=$$1} END {print s}'
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

    override lazy val rawInputs = Map("sc_test.do_prepare.input_file" -> createCannedFile("scatter",contents).toString,
    "sc_test.do_scatter.salt" -> salt)
  }

  object FileClobber extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      s"""task read_line {
        |  File in
        |  command { cat $${in} }
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
      "two.x.in" -> firstFile.toString,
      "two.y.in" -> secondFile.toString
    )
  }

  object FilePassingWorkflow extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      s"""task a {
        |  File in
        |  String out_name = "out"
        |
        |  command {
        |    cat $${in} > $${out_name}
        |  }
        |  RUNTIME
        |  output {
        |    File out = "out"
        |    File out_interpolation = "$${out_name}"
        |    String contents = read_string("$${out_name}")
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
      "file_passing.f" -> createCannedFile("canned", fileContents).toString
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
      s"""task a {
        |  File in
        |  String out_name = "out"
        |  String salt
        |
        |  command {
        |    # $${salt}
        |    echo "Something"
        |    cat $${in} > $${out_name}
        |  }
        |  RUNTIME
        |  output {
        |    File out = "out"
        |    File out_interpolation = "$${out_name}"
        |    String contents = read_string("$${out_name}")
        |    Array[String] stdoutContent = read_lines(stdout())
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
      "file_passing.f" -> createCannedFile("canned", fileContents).toString,
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
      s"""
        |task a {
        |  Array[String] array
        |  Map[String, String] map
        |
        |  command {
        |    echo $${sep = ' ' array} > concat
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
        "w.array_file" -> createCannedFile("array.txt", CannedArray).toString,
        "w.map_file" -> createCannedFile("map.txt", CannedMap).toString
      )
    }
  }

  object ArrayOfArrays extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      s"""task subtask {
        |  Array[File] a
        |  command {
        |    cat $${sep = " " a}
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
          WdlArray(WdlArrayType(WdlStringType), Seq(firstFile.toString, secondFile.toString).map(WdlString)),
          WdlArray(WdlArrayType(WdlStringType), Seq(thirdFile.toString, fourthFile.toString).map(WdlString))
        )
      )
    )
  }

  object CallCachingHashingWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      s"""task t {
        |  Int a
        |  Float b
        |  String c
        |  File d
        |
        |  command {
        |    echo "$${a}" > a
        |    echo "$${b}" > b
        |    echo "$${c}" > c
        |    cat $${d} > d
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
      "w.t.d" -> WdlFile(cannedFile.toString)
    )
  }

  object ExpressionsInInputs extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      s"""task echo {
        |  String inString
        |  command {
        |    echo $${inString}
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

  object WorkflowFailSlow extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      s"""
task shouldCompleteFast {
        |    Int a
        |    command {
        |        echo "The number was: $${a}"
        |    }
        |    output {
        |        Int echo = a
        |    }
        |}
        |
        |task shouldCompleteSlow {
        |    Int a
        |    command {
        |        echo "The number was: $${a}"
        |        # More than 1 so this should finish second
        |        sleep 2
        |    }
        |    output {
        |        Int echo = a
        |    }
        |}
        |
        |task failMeSlowly {
        |    Int a
        |    command {
        |        echo "The number was: $${a}"
        |        # Less than 2 so this should finish first
        |        sleep 1
        |        ./NOOOOOO
        |    }
        |    output {
        |        Int echo = a
        |    }
        |}
        |
        |task shouldNeverRun {
        |    Int a
        |    Int b
        |    command {
        |        echo "You can't fight in here - this is the war room $${a + b}"
        |    }
        |    output {
        |        Int echo = a
        |    }
        |}
        |
        |workflow wf {
        |    call shouldCompleteFast as A { input: a = 5 }
        |    call shouldCompleteFast as B { input: a = 5 }
        |
        |    call failMeSlowly as ohNOOOOOOOO { input: a = A.echo }
        |    call shouldCompleteSlow as C { input: a = B.echo }
        |
        |    call shouldNeverRun as D { input: a = ohNOOOOOOOO.echo, b = C.echo }
        |    call shouldCompleteSlow as E { input: a = C.echo }
        |}
      """.stripMargin

    val rawInputs = Map(
      "w.x.a" -> 5
    )
  }
}
