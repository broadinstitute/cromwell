package cromwell

import java.io.{File, FileWriter}
import java.util.UUID
import akka.pattern.ask
import akka.actor.ActorSystem
import akka.testkit.{TestFSMRef, TestActorRef, filterEvents}
import com.typesafe.config.ConfigFactory
import cromwell.binding._
import cromwell.binding.values.{WdlInteger, WdlValue}
import cromwell.engine.{WorkflowRunning, WorkflowSucceeded, WorkflowSubmitted, WorkflowActor}
import cromwell.engine.WorkflowActor._
import cromwell.engine.backend.local.LocalBackend

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps

object ThreeStepActorSpec {
  val Config =
    """
      |akka {
      |  loggers = ["akka.testkit.TestEventListener"]
      |  loglevel = "DEBUG"
      |  actor.debug.receive = on
      |}
    """.stripMargin

  val WdlSource =
    """
      |task ps {
      |  command {
      |    cat ${filename}
      |  }
      |  output {
      |    File procs = stdout()
      |  }
      |}
      |
      |task cgrep {
      |  command {
      |    grep '${pattern}' ${File in_file} | wc -l
      |  }
      |  output {
      |    Int count = read_int(stdout())
      |  }
      |}
      |
      |task wc {
      |  command {
      |    cat ${File in_file} | wc -l
      |  }
      |  output {
      |    Int count = read_int(stdout())
      |  }
      |}
      |
      |workflow three_step {
      |  call ps
      |  call cgrep {
      |    input: in_file=ps.procs
      |  }
      |  call wc {
      |    input: in_file=ps.procs
      |  }
      |}
    """.stripMargin

  val DummyProcessOutput =
    """
      |USER              PID  %CPU %MEM      VSZ    RSS   TT  STAT STARTED      TIME COMMAND
      |joeblaux         10344   4.5 20.6  7011548 3454616  ??  S    Mon06AM 275:26.10 /Applications/IntelliJ IDEA 14.app/Contents/MacOS/idea
      |joeblaux           883   2.2  0.5  2716336  85768   ??  S    Sun08AM   0:52.64 /Applications/Utilities/Terminal.app/Contents/MacOS/Terminal
      |_coreaudiod        577   1.9  0.1  2522428   9572   ??  Ss   Sun08AM  55:39.69 /usr/sbin/coreaudiod
      |_windowserver      130   1.6  1.9  4535136 319588   ??  Ss   Sun08AM 148:39.39 /System/Library/Frameworks/ApplicationServices.framework/Frameworks/CoreGraphics.framework/Resources/WindowServer -daemon
      |joeblaux          6440   1.4  2.2  4496164 362136   ??  S    Sun09PM  74:29.40 /Applications/Google Chrome.app/Contents/MacOS/Google Chrome
    """.stripMargin.trim

  val TestExecutionTimeout = 1500 milliseconds

  object Inputs {
    val Pattern = "three_step.cgrep.pattern"
    val DummyPsFileName = "three_step.ps.filename"
  }

  def createDummyPsFile: File = {
    val file = File.createTempFile("dummy_ps", ".out")
    val writer = new FileWriter(file)
    writer.write(DummyProcessOutput)
    writer.flush()
    writer.close()
    file
  }
}

class ThreeStepActorSpec extends CromwellSpec(ActorSystem("ThreeStepActorSpec", ConfigFactory.parseString(ThreeStepActorSpec.Config))) {
  private def buildWorkflowActor: TestActorRef[WorkflowActor] = {
    val (namespace, coercedInputs) = buildWorkflowBindingAndInputs()
    val props = WorkflowActor.props(UUID.randomUUID(), namespace, coercedInputs, new LocalBackend)
    TestActorRef(props, self, "ThreeStep")
  }

  private def buildWorkflowActorFsm() = {
    val (binding, coercedInputs) = buildWorkflowBindingAndInputs()
    TestFSMRef(new WorkflowActor(UUID.randomUUID(), binding, coercedInputs, new LocalBackend))
  }

  private def buildWorkflowBindingAndInputs() = {
    import ThreeStepActorSpec._
    val workflowInputs = Map(
      Inputs.Pattern ->"joeblaux",
      Inputs.DummyPsFileName -> createDummyPsFile.getAbsolutePath)

    val namespace = WdlNamespace.load(ThreeStepActorSpec.WdlSource)
    // This is a test and is okay with just throwing if coerceRawInputs returns a Failure.
    val coercedInputs = namespace.coerceRawInputs(workflowInputs).get
    val props = WorkflowActor.props(UUID.randomUUID(), namespace, coercedInputs, new LocalBackend)
    TestActorRef(props, self, "ThreeStep")
    (namespace, coercedInputs)
  }

  private def getCounts(outputs: Map[String, WdlValue], outputFqns: String*): Seq[Int] = {
    outputFqns.toSeq.map { outputFqn =>
      val wdlValue = outputs.getOrElse(outputFqn, throw new RuntimeException(s"Output $outputFqn not found"))
      wdlValue.asInstanceOf[WdlInteger].value.toInt
    }
  }

  import ThreeStepActorSpec._

  "A three step workflow" should {
      "best get to (three) steppin'" in {
      val fsm = buildWorkflowActorFsm()
      assert(fsm.stateName == WorkflowSubmitted)
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
}
