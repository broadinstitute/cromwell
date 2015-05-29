package cromwell.binding

import org.scalatest.{FlatSpec, Matchers}

class ThreeStepImportSpec extends FlatSpec with Matchers {
  val psTaskWdl = """
    |task ps {
    |  command {
    |    ps
    |  }
    |  output {
    |    File procs = stdout()
    |  }
    |}""".stripMargin

  val cgrepTaskWdl = """
    |task cgrep {
    |  command {
    |    grep '${pattern}' ${File in_file} | wc -l
    |  }
    |  output {
    |    Int count = read_int(stdout())
    |  }
    |}""".stripMargin

  val wcTaskWdl = """
    |task wc {
    |  command {
    |    cat ${File in_file} | wc -l
    |  }
    |  output {
    |    Int count = read_int(stdout())
    |  }
    |}""".stripMargin

  val workflowWdl = """
    |import "ps"
    |import "cgrep"
    |import "wc"
    |
    |workflow three_step {
    |  call ps
    |  call cgrep {
    |    input: in_file=ps.procs
    |  }
    |  call wc {
    |    input: in_file=ps.procs
    |  }
    |}""".stripMargin

  def resolver(importUri: String): String = {
    importUri match {
      case "ps" => psTaskWdl
      case "cgrep" => cgrepTaskWdl
      case "wc" => wcTaskWdl
      case _ => throw new RuntimeException(s"Can't resolve $importUri")
    }
  }

  val binding = WdlBinding.process(workflowWdl, resolver _)

  "WDL file with imports" should "Have 4 executables (3 tasks, 1 workflow)" in {
    binding.executables.size shouldEqual 4
  }
  it should "Have 3 tasks" in {
    binding.tasks.size shouldEqual 3
  }
  it should "Have 3 task ASTs" in {
    binding.taskAsts.size shouldEqual 3
  }
  it should "Have 0 local tasks" in {
    binding.localTasks.size shouldEqual 0
  }
  it should "Have 0 local task ASTs" in {
    binding.localTaskAsts.size shouldEqual 0
  }
  it should "Have 3 imported tasks" in {
    binding.importedTasks.size shouldEqual 3
  }
  it should "Have 3 imported task ASTs" in {
    binding.importedTaskAsts.size shouldEqual 3
  }
  it should "Have 1 workflow" in {
    binding.workflows.size shouldEqual 1
  }
  it should "Have 1 local workflow" in {
    binding.localWorkflows.size shouldEqual 1
  }
  it should "Have 0 imported workflow" in {
    binding.importedWorkflows.size shouldEqual 0
  }
  it should "Have 3 imported WdlBindings" in {
    binding.importedBindings.size shouldEqual 3
  }
  it should "Have tasks with the names 'ps', 'cgrep' and 'wc'" in {
    binding.tasks map {_.name} shouldEqual Seq("ps", "cgrep", "wc")
  }
  it should "Throw an exception if the import resolver fails to resolve an import" in {
    def badResolver(s: String): String = {
      throw new RuntimeException(s"Can't Resolve")
    }
    try {
      val badBinding = WdlBinding.process(workflowWdl, badResolver _)
      fail("Expecting an exception to be thrown when using badResolver")
    } catch {
      case _: RuntimeException =>
    }
  }
}
