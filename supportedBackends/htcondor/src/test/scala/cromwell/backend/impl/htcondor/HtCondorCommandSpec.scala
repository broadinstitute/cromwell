package cromwell.backend.impl.htcondor

import better.files._

import org.scalatest.{Matchers, WordSpecLike}

class HtCondorCommandSpec extends WordSpecLike with Matchers{
  val attributes = Map("executable" -> "test.sh", "input" -> "/temp/test", "error"->"stderr")
  val resultAttributes = List("executable=test.sh","input=/temp/test","error=stderr", "queue")
  val htCondorCommands = new HtCondorCommands

  "submitCommand method" should {
    "return submit file with content passed to it" in {
      val dir = File.newTemp()
      val command = htCondorCommands.generateSubmitFile(dir.path,attributes)
      val file = dir
      resultAttributes shouldEqual  dir.lines.toList
      dir.delete()
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