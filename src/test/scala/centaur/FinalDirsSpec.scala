package centaur

import java.nio.file.{Paths, Path}

import cats.data.Validated.{Invalid, Valid}
import centaur.test.formulas.TestFormulas
import centaur.test.workflow.Workflow
import org.scalatest.{FlatSpec, Matchers, ParallelTestExecution}

object FinalDirsSpec {
  val FinalDir = Paths.get("src/main/resources/finalCopy")
  val OutputsDirTest = FinalDir.resolve("final_workflow_outputs.test")
  val LogDirTest = FinalDir.resolve("final_workflow_log_dir.test")
  val CallLogsDirTest = FinalDir.resolve("final_call_logs_dir.test")
}

class FinalDirsSpec extends FlatSpec with Matchers with ParallelTestExecution {
  import FinalDirsSpec._

  "finalOutputs" should "place files in output dir" in testFinalOutputs(OutputsDirTest, "final_workflow_outputs_dir")

  "finalDirs" should "place log files in output dir when requested" in testFinalOutputs(LogDirTest, "final_workflow_log_dir")

  "finalOutputs" should "place call files in output dir when requested" in testFinalOutputs(CallLogsDirTest, "final_call_logs_dir")

  def testFinalOutputs(path: Path, option: String) = 
    Workflow.fromPath(path) match {
      case Valid(w) => TestFormulas.runFinalDirsWorkflow(w, option).run.get
      case Invalid(e) => fail(s"Could not read logs test:\n -${e.toList.mkString("\n-")}")
    }

}
