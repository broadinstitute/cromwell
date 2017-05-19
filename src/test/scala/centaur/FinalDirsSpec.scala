package centaur
import centaur.test.{CheckFiles, JesCheckFiles, LocalCheckFiles, TestOptions, Test}
import java.nio.file.{Path, Paths}
import java.io.File
import cats._, cats.implicits._

import cats.data.Validated.{Invalid, Valid}
import centaur.api.CentaurCromwellClient
import centaur.test.formulas.TestFormulas
import centaur.test.standard.StandardTestCase
import centaur.test.standard.StandardTestFormat.WorkflowSuccessTest
import centaur.test.workflow.Workflow
import org.scalatest.{FlatSpec, Matchers, ParallelTestExecution}
import com.google.cloud.storage._

object FinalDirsSpec {
  val FinalDir = Paths.get("src/main/resources/finalCopy")

  val OutputsDirLocal= FinalDir.resolve("final_workflow_outputs_local.test")
  val OutputsDirJes= FinalDir.resolve("final_workflow_outputs_jes.test")
  val OutputsDirTes= FinalDir.resolve("final_workflow_outputs_local.test")

  val LogDirLocal = FinalDir.resolve("final_workflow_log_dir_local.test")
  val LogDirJes = FinalDir.resolve("final_workflow_log_dir_jes.test")
  val LogDirTes = FinalDir.resolve("final_workflow_log_dir_local.test")

  val CallLogsDirLocal = FinalDir.resolve("final_call_logs_dir_local.test")
  val CallLogsDirJes = FinalDir.resolve("final_call_logs_dir_jes.test")
  val CallLogsDirTes = FinalDir.resolve("final_call_logs_dir_tes.test")
}


class FinalDirsSpec extends FlatSpec with Matchers with ParallelTestExecution {
  import FinalDirsSpec._

  case class Backend(id: String, checkFiles: CheckFiles)
  val Jes = Backend("jes", JesCheckFiles(StorageOptions.getDefaultInstance.getService))
  val Local = Backend("local", LocalCheckFiles())
  val Tes = Backend("tes", LocalCheckFiles())

  val cromwellBackends = CentaurCromwellClient.backends.get.supportedBackends.map(_.toLowerCase)

  "final Outputs on local backend" should "place files in output dir" in
    testFinalOutputs(OutputsDirLocal, "final_workflow_outputs_dir", Local)

  "final log dirs on local backend" should "place log files in output dir when requested" in
    testFinalOutputs(LogDirLocal, "final_workflow_log_dir", Local)

  "final call logs dir in local" should "place call files in output dir when requested" in 
    testFinalOutputs(CallLogsDirLocal, "final_call_logs_dir", Local)


  "final Logs dir on jes backend" should "place log files in output dir when requested" in
    testFinalOutputs(LogDirJes, "final_workflow_log_dir", Jes)

  "final outputs on jes backend" should "place files in output dir" in 
    testFinalOutputs(OutputsDirJes, "final_workflow_outputs_dir", Jes)

  "final Call logs dir in Jes" should "place call files in GCS dir when requested" in 
    testFinalOutputs(CallLogsDirJes, "final_call_logs_dir", Jes)
    

  "final Call logs dir in Tes" should "place call files in local dir when requested" in 
    testFinalOutputs(OutputsDirTes, "final_workflow_outputs_dir", Tes)

  "final Logs dir on Tes backend" should "place log files in output dir when requested" in
    testFinalOutputs(LogDirTes, "final_workflow_log_dir", Tes)

  "final Call logs dir in Tes" should "place call files in local tes logs dir when requested" in 
    testFinalOutputs(CallLogsDirTes, "final_call_logs_dir", Tes)


  def testFinalOutputs(path: Path, option: String, backend: Backend) =
    (Workflow.fromPath(path), cromwellBackends.contains(backend.id)) match {

      case (Valid(w), true) =>
        TestFormulas.runFinalDirsWorkflow(w, option, backend.checkFiles).run.get

      case (Invalid(e), _) => fail(s"Could not read logs test:\n -${e.toList.mkString("\n-")}")

      case (Valid(_), _) => println(s"didn't find ${backend.id} in $cromwellBackends"); ignore
    }
}
