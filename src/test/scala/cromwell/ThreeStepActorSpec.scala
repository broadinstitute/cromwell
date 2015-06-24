package cromwell

import java.io.{File, FileWriter}

import akka.pattern.ask
import akka.testkit.TestFSMRef
import cromwell.CromwellSpec.DockerTest
import cromwell.binding._
import cromwell.binding.values.{WdlInteger, WdlValue}
import cromwell.engine.workflow.WorkflowActor
import WorkflowActor._
import cromwell.engine._
import cromwell.engine.backend.local.LocalBackend
import cromwell.util.SampleWdl

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps


object ThreeStepActorSpec {
  val DummyProcessOutput =
    """
      |USER              PID  %CPU %MEM      VSZ    RSS   TT  STAT STARTED      TIME COMMAND
      |joeblaux         10344   4.5 20.6  7011548 3454616  ??  S    Mon06AM 275:26.10 /Applications/IntelliJ IDEA 14.app/Contents/MacOS/idea
      |joeblaux           883   2.2  0.5  2716336  85768   ??  S    Sun08AM   0:52.64 /Applications/Utilities/Terminal.app/Contents/MacOS/Terminal
      |_coreaudiod        577   1.9  0.1  2522428   9572   ??  Ss   Sun08AM  55:39.69 /usr/sbin/coreaudiod
      |_windowserver      130   1.6  1.9  4535136 319588   ??  Ss   Sun08AM 148:39.39 /System/Library/Frameworks/ApplicationServices.framework/Frameworks/CoreGraphics.framework/Resources/WindowServer -daemon
      |joeblaux          6440   1.4  2.2  4496164 362136   ??  S    Sun09PM  74:29.40 /Applications/Google Chrome.app/Contents/MacOS/Google Chrome
    """.stripMargin.trim

  val TestExecutionTimeout = 5000 milliseconds

  def createDummyPsFile: File = {
    val file = File.createTempFile("dummy_ps", ".out")
    val writer = new FileWriter(file)
    writer.write(DummyProcessOutput)
    writer.flush()
    writer.close()
    file
  }
}

class ThreeStepActorSpec extends CromwellTestkitSpec("ThreeStepActorSpec") {
  import ThreeStepActorSpec._

  private def buildWorkflowActorFsm(runtime: String) = {
    import SampleWdl.FauxThreeStep.InputKeys._
    val workflowInputs = Map(
      Pattern -> "joeblaux",
      DummyPsFile -> createDummyPsFile.getAbsolutePath,
      DummyPs2File -> createDummyPsFile.getAbsolutePath,
      DummyPs3File -> createDummyPsFile.getAbsolutePath)

    val namespace = WdlNamespace.load(SampleWdl.FauxThreeStep.wdlSource(runtime))
    // This is a test and is okay with just throwing if coerceRawInputs returns a Failure.
    val coercedInputs = namespace.coerceRawInputs(workflowInputs).get
    val descriptor = WorkflowDescriptor(namespace, coercedInputs)
    TestFSMRef(new WorkflowActor(descriptor, new LocalBackend))
  }

  private def getCounts(outputs: Map[String, WdlValue], outputFqns: String*): Seq[Int] = {
    outputFqns.toSeq.map { outputFqn =>
      val wdlValue = outputs.getOrElse(outputFqn, throw new RuntimeException(s"Output $outputFqn not found"))
      wdlValue.asInstanceOf[WdlInteger].value.toInt
    }
  }

  private def runAndAssertCorrectness(runtime: String = ""): Unit = {
    val fsm = buildWorkflowActorFsm(runtime)
    assert(fsm.stateName == WorkflowSubmitted)
    startingCallsFilter("cgrep", "wc") {
      fsm ! Start
      within(TestExecutionTimeout) {
        awaitCond(fsm.stateName == WorkflowRunning)
        awaitCond(fsm.stateName == WorkflowSucceeded)
        val outputs = Await.result(fsm.ask(GetOutputs).mapTo[WorkflowOutputs], 5 seconds)
        val Seq(cgrepCount, wcCount) = getCounts(outputs, "three_step.cgrep.count", "three_step.wc.count")
        cgrepCount shouldEqual 3
        wcCount shouldEqual 6
      }
    }
  }

  "A three step workflow" should {
    "best get to (three) steppin'" in {
      runAndAssertCorrectness()
    }
  }

  "A Dockerized three step workflow" should {
    "best get to Dockerized (three) steppin'" taggedAs DockerTest in {
      runAndAssertCorrectness(
        """
          |runtime {
          |  docker: "ubuntu:latest"
          |}
        """.stripMargin)
    }
  }
}
