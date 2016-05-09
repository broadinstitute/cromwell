package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import wdl4s.values.WdlString
import cromwell.util.SampleWdl

import scala.language.postfixOps

class GlobbingWorkflowSpec extends CromwellTestkitSpec {
  def doTheTest(runtime: String = "") = {
    val outputs = runWdl(
      sampleWdl = SampleWdl.GlobtasticWorkflow,
      eventFilter = EventFilter.info(pattern = s"starting calls: w.B", occurrences = 1),
      runtime = runtime
    )

    // The order in which files glob is apparently not guaranteed, so accept any permutation.
    val permutations = for {
      permutation <- Seq('a', 'b', 'c').permutations
    } yield WdlString(permutation.mkString("\n"))

    val actual = outputs.get("w.B.B_out").get.wdlValue
    permutations collectFirst { case s: WdlString if s == actual => s } should not be empty
  }

  val newline = System.lineSeparator
  "A workflow with globbed outputs" should {
    "run properly" ignore doTheTest()
    "run properly in a Docker environment" taggedAs DockerTest ignore doTheTest("""runtime { docker: "ubuntu:latest" }""")
  }
}
