package cromwell.backend.impl.htcondor

import better.files._

import org.scalatest.{Matchers, WordSpecLike}

class HtCondorCommandSpec extends WordSpecLike with Matchers {
  private val attributes = Map("executable" -> "test.sh", "input" -> "/temp/test", "error"->"stderr")
  private val resultAttributes = List("executable=test.sh","input=/temp/test","error=stderr", "spec1", "spec2", "queue")
  private val htCondorCommands = new HtCondorCommands
  private val nativeSpecs = Option(Array("spec1", "spec2"))

  "submitCommand method" should {
    "return submit file with content passed to it" in {
      val file = File.newTemporaryFile()
      val command = htCondorCommands.generateSubmitFile(file.path, attributes, nativeSpecs)
      resultAttributes shouldEqual file.lines.toList
      file.delete()
      command shouldEqual s"condor_submit ${file.path}"
    }
  }

  "statusCommand method" should {
    "return status command" in {
      val command = HtCondorCommands.generateJobStatusCommand("96.0")
      command shouldEqual s"condor_q 96.0 -autoformat JobStatus"
    }
  }
}