package wdl4s

import better.files.File
import org.scalatest.{FlatSpec, Matchers}


class WorkflowImportsSpec extends FlatSpec with Matchers {

  def addAndGetFile(name: String, source: String): String = {
    val tempFile = File.newTemporaryFile(s"${name}", ".wdl", Option(wdlDirectory)) write source
    tempFile.name
  }

  private val wdlDirectory = File.newTemporaryDirectory("imports_dir")

  val echoHelloWdl =
    """
      |task inItalian {
      |  String x = "ciao"
      |  command {
      |    echo "${x}"
      |  }
      |}
      |
      |task inSpanish {
      |  String x = "hola"
      |  command {
      |    echo "${x}"
      |  }
      |}
      |
      |task inFrench {
      |  String x = "bonjour"
      |  command {
      |    echo "${x}"
      |  }
      |}
      |
      |workflow echoHello {}
    """.stripMargin

  val echoHelloWdlFile = addAndGetFile("echoHello", echoHelloWdl)

  val basicWdl =
    s"""
      |import "$echoHelloWdlFile"
      """.stripMargin +
    """
      |
      |task ls {
      |  command {
      |    ls -l
      |  }
      |  output {
      |    String fileList = read_string(stdout())
      |  }
      |}
      |
      |task pwd {
      |   command {
      |     pwd
      |   }
      |   output {
      |     String current = read_string(stdout())
      |   }
      |}
      |
      |workflow basic {
      |  call ls
      |  call pwd
      |  output {
      |       ls.fileList
      |   }
      |}
    """.stripMargin

  val printNumsWdl =
    s"""
      |import "$echoHelloWdlFile" as multilingualEcho
      |
      """.stripMargin +
    """
      |task ls {
      |  command {
      |    ls
      |  }
      |}
      |
      |
      |task print1 {
      |  Int x = 10
      |  command {
      |    for i in `seq 1 ${x}`
      |      do
      |        echo $i
      |      done
      |  }
      |}
      |
      |task print2 {
      |  Int x = 20
      |  command {
      |    for i in `seq 1 ${x}`
      |      do
      |        echo $i
      |      done
      |  }
      |}
      |
      |task print3 {
      |  Int x = 30
      |  command {
      |    for i in `seq 1 ${x}`
      |      do
      |        echo $i
      |      done
      |  }
      |}
    """.stripMargin

  val snoozeWdl =
    """
      |task sleep {
      |  command { sleep 1 }
      |
      |}
      |task sleep2 {
      |  command { sleep 2 }
      |}
      |
      |task sleep3 {
      |  command { sleep 3 }
      |}
    """.stripMargin

  val threeStepWdl =
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
      |  command {
      |    grep '${pattern}' ${in_file} | wc -l
      |  }
      |  output {
      |    Int count = read_int(stdout())
      |  }
      |  runtime {docker: "ubuntu:latest"}
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
      |  runtime {docker: "ubuntu:latest"}
      |}
      |
      |workflow three_step {
      |  call ps
      |  call cgrep { input: in_file=ps.procs }
      |  call wc { input: in_file=ps.procs }
      |}
    """.stripMargin


  val threeStepWdlWithImports =
    s"""
      |import "$echoHelloWdlFile" as funEcho
      |
      """.stripMargin + threeStepWdl

  val basicWdlImportFile = addAndGetFile("basic", basicWdl)
  val threeStepWdlImportFile = addAndGetFile("threestep", threeStepWdlWithImports)
  val snoozeWdlImportFile = addAndGetFile("snooze", snoozeWdl)
  val printNumsWdlImportFile = addAndGetFile("printNums", printNumsWdl)

  def noExtension(fileName: String) = fileName.replace(".wdl","")


  val imports =
    s"""
       |import "$basicWdlImportFile"
       |import "$threeStepWdlImportFile" as classicThreeStep
       |import "$snoozeWdlImportFile" as trySleep
       |import "$printNumsWdlImportFile"
       |
    """.stripMargin

  val primaryWorkflow =
    s"""
    |
    |task testCaseTask {
    |  command {
    |     echo "Ruchi's birthday: 01/19"
    |  }
    |}
    |
    |workflow testCases {
    |  call ${noExtension(basicWdlImportFile)}.ls as ls1
    |  call ${noExtension(printNumsWdlImportFile)}.ls as ls2
    |
    |  #call ${noExtension(basicWdlImportFile)}.pwd as soBasic
    |
    |  #call ${noExtension(printNumsWdlImportFile)}.print1
    |  #call ${noExtension(printNumsWdlImportFile)}.print2 as printingFun
    |
    |
    |  call classicThreeStep.ps
    |  call classicThreeStep.ps as psAgain
    |
    |  call trySleep.sleep
    |  call trySleep.sleep2 as sleepMore
    |
    |  call ${noExtension(basicWdlImportFile)}.${noExtension(echoHelloWdlFile)}.inFrench
    |  call classicThreeStep.funEcho.inSpanish
    |  call classicThreeStep.funEcho.inSpanish as inPortugese
    |  call ${noExtension(printNumsWdlImportFile)}.multilingualEcho.inItalian
    |
    |}
    |""".stripMargin

  val wdlWithImports = imports + primaryWorkflow

  val namespace = WdlNamespaceWithWorkflow.load(wdlWithImports, wdlDirectory)

  "WDL file with imports" should "Have 1 task (remaining tasks are in separate namespaces)" in {
    namespace.tasks.size shouldEqual 1
  }

  it should "Have 4 imported WdlNamespaces" in {
    namespace.namespaces.size shouldEqual 4
  }

  it should "import a WDL file (with alias) and be able to reference its tasks by FQN" in {
    namespace.resolve("classicThreeStep.ps").size shouldEqual 1
  }
  it should "import a WDL file (with no alias) and be able to reference its tasks by FQN" in {
    namespace.resolve(s"${noExtension(printNumsWdlImportFile)}.print1").size shouldEqual 1
  }
  it should "import two WDL file (with clashing task names) and be able to reference all tasks by FQN" in {
    val clashingTaskNames = Seq(namespace.resolve(s"${noExtension(basicWdlImportFile)}.ls"),
                            namespace.resolve(s"${noExtension(printNumsWdlImportFile)}.ls"))

    clashingTaskNames.size shouldEqual 2
  }

  def deleteTempFiles() = wdlDirectory.delete(swallowIOExceptions = true)
  deleteTempFiles()
}
