package cromwell

import akka.testkit.EventFilter
import cromwell.util.SampleWdl
import wdl4s._
import wdl4s.values.WdlFloat

class FileSizeWorkflowSpec extends CromwellTestkitSpec {

  object FileSize extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """task file_size {
        |  File input_file
        |
        |  command {
        |    echo "this file is 22 bytes" > created_file
        |  }
        |
        |  output {
        |    Float input_file_size = size(input_file)
        |    Float created_file_size = size("created_file")
        |    Float created_file_size_in_k = size("created_file", "K")
        |    Float created_file_size_in_kb = size("created_file", "KB")
        |    Float created_file_size_in_m = size("created_file", "M")
        |    Float created_file_size_in_mb = size("created_file", "MB")
        |    Float created_file_size_in_g = size("created_file", "G")
        |    Float created_file_size_in_gb = size("created_file", "GB")
        |    Float created_file_size_in_t = size("created_file", "T")
        |    Float created_file_size_in_tb = size("created_file", "TB")
        |    Float created_file_size_in_ki = size("created_file", "Ki")
        |    Float created_file_size_in_kib = size("created_file", "KiB")
        |    Float created_file_size_in_mi = size("created_file", "Mi")
        |    Float created_file_size_in_mib = size("created_file", "MiB")
        |    Float created_file_size_in_gi = size("created_file", "Gi")
        |    Float created_file_size_in_gib = size("created_file", "GiB")
        |    Float created_file_size_in_ti = size("created_file", "Ti")
        |    Float created_file_size_in_tib = size("created_file", "TiB")
        |  }
        |}
        |
        |workflow wf {
        |  call file_size
        |}
      """.stripMargin
    override val rawInputs: WorkflowRawInputs = Map(
      "wf.file_size.input_file" -> createCannedFile("canned", "arbitrary content").getAbsolutePath
    )
  }

  "size engine function" should {
    "return the size of a file in bytes or convert it to the specified unit" in {
      val outputs = Map(
        "wf.file_size.input_file_size" -> WdlFloat(17),
        "wf.file_size.created_file_size" -> WdlFloat(22),
        "wf.file_size.created_file_size_in_k" -> WdlFloat(0.022),
        "wf.file_size.created_file_size_in_kb" -> WdlFloat(0.022),
        "wf.file_size.created_file_size_in_m" -> WdlFloat(0.000022),
        "wf.file_size.created_file_size_in_mb" -> WdlFloat(0.000022),
        "wf.file_size.created_file_size_in_g" -> WdlFloat(0.000000022),
        "wf.file_size.created_file_size_in_gb" -> WdlFloat(0.000000022),
        "wf.file_size.created_file_size_in_t" -> WdlFloat(0.000000000022),
        "wf.file_size.created_file_size_in_tb" -> WdlFloat(0.000000000022),
        "wf.file_size.created_file_size_in_ki" -> WdlFloat(0.021484375),
        "wf.file_size.created_file_size_in_kib" -> WdlFloat(0.021484375),
        "wf.file_size.created_file_size_in_mi" -> WdlFloat(2.09808349609375E-5),
        "wf.file_size.created_file_size_in_mib" -> WdlFloat(2.09808349609375E-5),
        "wf.file_size.created_file_size_in_gi" -> WdlFloat(2.0489096641540527E-8),
        "wf.file_size.created_file_size_in_gib" -> WdlFloat(2.0489096641540527E-8),
        "wf.file_size.created_file_size_in_ti" -> WdlFloat(2.000888343900442E-11),
        "wf.file_size.created_file_size_in_tib" -> WdlFloat(2.000888343900442E-11)
      )
      runWdlAndAssertOutputs(
        sampleWdl = FileSize,
        eventFilter = EventFilter.info(pattern = s"transition from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1),
        expectedOutputs = outputs
      )
    }
  }

}
