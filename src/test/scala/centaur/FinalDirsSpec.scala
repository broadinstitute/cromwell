package centaur

import java.nio.file.{Paths, Path}

import java.io.File
import cats.data.Validated.{Invalid, Valid}
import centaur.test.formulas.TestFormulas
import centaur.test.workflow.Workflow
import org.scalatest.{FlatSpec, Matchers, ParallelTestExecution}

object FinalDirsSpec {
  val FinalDir = Paths.get("src/main/resources/finalCopy")
  val OutputsDirTest = FinalDir.resolve("final_workflow_outputs_local.test")
  val LogDirTest = FinalDir.resolve("final_workflow_log_dir_local.test")
  val CallLogsDirTest = FinalDir.resolve("final_call_logs_dir_local.test")
}

class FinalDirsSpec extends FlatSpec with Matchers with ParallelTestExecution {
  import FinalDirsSpec._

  "finalOutputs" should "place files in output dir" in testFinalOutputs(OutputsDirTest, "final_workflow_outputs_dir")

  "finalDirs" should "place log files in output dir when requested" in testFinalOutputs(LogDirTest, "final_workflow_log_dir")

  "finalOutputs" should "place call files in output dir when requested" in testFinalOutputs(CallLogsDirTest, "final_call_logs_dir")

  def testFinalOutputs(path: Path, option: String) = 
    Workflow.fromPath(path) match {
      case Valid(w) => TestFormulas.runFinalDirsWorkflow(w, option, checkDirectorySize, deleteLocalFiles(_)).run.get
      case Invalid(e) => fail(s"Could not read logs test:\n -${e.toList.mkString("\n-")}")
    }

  def checkDirectorySize: String => Int = {s =>
    (new java.io.File(s)).listFiles.size
  }


  def deleteLocalFiles(directory: String) {
    val dir = new java.io.File(directory)

    def deleteDir(d: File) {
      d.listFiles.toList.foreach { f =>

        //directories must be empty to be deleted
        if (f.isDirectory) 
          deleteDir(f)

        f.delete
      }
    }

    deleteDir(dir)
  }
}
