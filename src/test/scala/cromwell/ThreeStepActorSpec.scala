package cromwell

import java.io.{File, FileWriter}
import java.util.UUID

import akka.actor.ActorSystem
import akka.testkit.filterEvents
import akka.testkit.TestActorRef
import com.typesafe.config.ConfigFactory
import cromwell.binding._
import cromwell.binding.values.{WdlInteger, WdlString}
import cromwell.engine.WorkflowActor._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.{SymbolStore, WorkflowActor}

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
      |    File procs = "stdout"
      |  }
      |}
      |
      |task cgrep {
      |  command {
      |    grep '${pattern}' ${File in_file} | wc -l
      |  }
      |  output {
      |    Int count = read_int("stdout")
      |  }
      |}
      |
      |task wc {
      |  command {
      |    cat ${File in_file} | wc -l
      |  }
      |  output {
      |    Int count = read_int("stdout")
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

  val TestExecutionTimeout = 500 milliseconds

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

  val binding = WdlBinding.process(ThreeStepActorSpec.WdlSource)

  private def buildWorkflowActor: TestActorRef[WorkflowActor] = {
    import ThreeStepActorSpec._
    val workflowInputs = Map(
      Inputs.Pattern -> WdlString("joeblaux"),
      // TODO would this work as a WdlFile?
      Inputs.DummyPsFileName -> WdlString(createDummyPsFile.getAbsolutePath))

    val binding = WdlBinding.process(ThreeStepActorSpec.WdlSource)
    val props = WorkflowActor.buildWorkflowActorProps(UUID.randomUUID(), binding, workflowInputs)
    TestActorRef(props, "ThreeStep")
  }

  private def getCounts(symbolStore: SymbolStore, outputFqns: String*): Seq[Int] = {
    outputFqns.toSeq.map { outputFqn =>
      val maybeOutput = symbolStore.getOutputByFullyQualifiedName(outputFqn)
      val symbolStoreEntry = maybeOutput.getOrElse(throw new RuntimeException("No symbol store entry found!"))
      val wdlValue = symbolStoreEntry.wdlValue.getOrElse(throw new RuntimeException("No workflow output found!"))
      wdlValue.asInstanceOf[WdlInteger].value.toInt
    }
  }

  import ThreeStepActorSpec._

//  "A three step workflow" should {
//    "best get to (three) steppin'" in {
//      within(TestExecutionTimeout) {
//        filterEvents(startingCallsFilter("ps"), startingCallsFilter("cgrep", "wc")) {
//          buildWorkflowActor ! Start(new LocalBackend)
//          expectMsgPF() {
//            case Started => ()
//          }
//          expectMsgPF() {
//            case Failed(t) =>
//              fail(t)
//            case Done(symbolStore) =>
//              val Seq(cgrepCount, wcCount) = getCounts(symbolStore, "three_step.cgrep.count", "three_step.wc.count")
//              cgrepCount shouldEqual 3
//              wcCount shouldEqual 6
//          }
//        }
//      }
//    }
//  }
}
