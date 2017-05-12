package centaur
import centaur.test.CheckFiles

import java.nio.file.{Paths, Path}

import java.io.File
import centaur.test.{LocalCheckFiles, JesCheckFiles}
import cats.data.Validated.{Invalid, Valid}
import centaur.test.formulas.TestFormulas
import centaur.test.workflow.Workflow
import org.scalatest.{FlatSpec, Matchers, ParallelTestExecution}
import com.google.cloud.storage._

object FinalDirsSpec {
  val FinalDir = Paths.get("src/main/resources/finalCopy")

  val OutputsDirLocalTest = FinalDir.resolve("final_workflow_outputs_local.test")
  val OutputsDirJesTest = FinalDir.resolve("final_workflow_outputs_jes.test")

  val LogDirLocalTest = FinalDir.resolve("final_workflow_log_dir_local.test")
  val LogDirJesTest = FinalDir.resolve("final_workflow_log_dir_jes.test")

  val CallLogsDirLocalTest = FinalDir.resolve("final_call_logs_dir_local.test")
  val CallLogsDirJesTest = FinalDir.resolve("final_call_logs_dir_jes.test")
}


class FinalDirsSpec extends FlatSpec with Matchers with ParallelTestExecution {
  import FinalDirsSpec._

  "final Outputs on local backend" should "place files in output dir" in 
    testFinalOutputs(OutputsDirLocalTest, "final_workflow_outputs_dir", LocalCheckFiles())

  "final log dirs on local backend" should "place log files in output dir when requested" in 
    testFinalOutputs(LogDirLocalTest, "final_workflow_log_dir", LocalCheckFiles())

  "final call logs dir in local" should "place call files in output dir when requested" in 
    testFinalOutputs(CallLogsDirLocalTest, "final_call_logs_dir", LocalCheckFiles())

  "final Logs dir on jes backend" should "place log files in output dir when requested" in 
    testFinalOutputs(LogDirJesTest, "final_workflow_log_dir", JesCheckFiles(StorageOptions.getDefaultInstance.getService))

  "final outputs on jes backend" should "place files in output dir" in 
    testFinalOutputs(OutputsDirJesTest, "final_workflow_outputs_dir", JesCheckFiles(StorageOptions.getDefaultInstance.getService))

  "final Call logs dir in Jes" should "place call files in GCS dir when requested" in 
    testFinalOutputs(CallLogsDirJesTest, "final_call_logs_dir", JesCheckFiles(StorageOptions.getDefaultInstance.getService))

  def testFinalOutputs(path: Path, option: String, backend: CheckFiles) = 
    Workflow.fromPath(path) match {
      case Valid(w) => TestFormulas.runFinalDirsWorkflow(w, option, backend).run.get
      case Invalid(e) => fail(s"Could not read logs test:\n -${e.toList.mkString("\n-")}")
    }


}
