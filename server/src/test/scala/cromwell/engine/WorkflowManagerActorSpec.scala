
package cromwell.engine

import akka.actor.ActorSystem
import com.typesafe.config.ConfigValueFactory
import cromwell.core._
import cromwell.engine.workflow.WorkflowDescriptorBuilderForSpecs
import cromwell.util.{SampleWdl, WorkflowImport}
import cromwell.{CromwellTestKitSpec, CromwellTestKitWordSpec}
import org.scalatest.enablers.Retrying
import wom.core.{ExecutableInputMap, WorkflowSource}

import scala.concurrent.duration._
import scala.language.postfixOps

class WorkflowManagerActorSpec extends CromwellTestKitWordSpec with WorkflowDescriptorBuilderForSpecs {
  override implicit val actorSystem: ActorSystem = system

  "A WorkflowManagerActor" should {

    // this test lol
    "run workflows in the correct directory" in {
      val outputs = runWdl(sampleWdl = SampleWdl.CurrentDirectory, testActorName = "TestCromwellRootActor-whereami")

      val outputName = "wf_whereami.whereami.pwd"
      val salutation = outputs(outputName)
      val actualOutput = salutation.valueString.trim
      actualOutput should endWith("/call-whereami/execution")
    }

    "pick up workflows" in {
      case class SubWorkflows(naptime: FiniteDuration) extends SampleWdl {
        override def workflowSource(runtime: String): WorkflowSource = root

        override val rawInputs: ExecutableInputMap = Map("parent.naptime" -> naptime.toSeconds.toInt)

        val root: WorkflowSource =
          """
            |version 1.0
            |
            |import "subsub.wdl" as sub
            |
            |workflow parent {
            |  input {
            |    Int naptime
            |  }
            |  scatter (i in range(3)) {
            |    call sub.sub { input: naptime = naptime }
            |  }
            |}
          """.stripMargin.trim

        val sub: WorkflowSource =
          """
            |version 1.0
            |
            |task snooze {
            |  input {
            |    Int naptime
            |  }
            |  command {
            |    echo "zzzz"; sleep ~{naptime}
            |  }
            |}
            |
            |workflow sub {
            |  input {
            |    Int naptime
            |  }
            |  call snooze { input: naptime = naptime }
            |}
            |
        """.stripMargin.trim

        override val imports: Option[Set[WorkflowImport]] =
          Option(
            Set(
              WorkflowImport(name = "subsub.wdl", content = sub)
            )
          )
      }

      val config = CromwellTestKitSpec.NooPServiceActorConfig.
        withValue("system.max-concurrent-workflows", ConfigValueFactory.fromAnyRef(2)).
        withValue("system.new-workflow-poll-rate", ConfigValueFactory.fromAnyRef(1))

      val rootActor = buildCromwellRootActor(config = config, actorName = "TestCromwellRootActor-pickup")
      val serviceRegistryActor = rootActor.underlyingActor.serviceRegistryActor

      val firstSources = SubWorkflows(naptime = 60 seconds).asWorkflowSources()

      def waitForState(workflowId: WorkflowId, state: WorkflowState): Unit = {
        eventually { verifyWorkflowState(serviceRegistryActor, workflowId, state) } (config = defaultPatience, pos = implicitly[org.scalactic.source.Position], retrying = implicitly[Retrying[Unit]])
      }

      val firstWorkflowId = rootActor.underlyingActor.submitWorkflow(firstSources)
      waitForState(firstWorkflowId, WorkflowRunning)
      // Give the subworkflows plenty of time to start running.
      Thread.sleep(10 * 1000L)

      // Submit a second workflow.
      val secondSources = SubWorkflows(naptime = 0 seconds).asWorkflowSources()
      val secondWorkflowId = rootActor.underlyingActor.submitWorkflow(secondSources)
      // Wait a long time so that there would have been many root workflow pickup sweeps.
      Thread.sleep(30 * 1000L)
      // Despite the long wait this workflow should not be running because the first workflow has more than
      // max-concurrent-workflows root and sub workflows running.
      waitForState(secondWorkflowId, WorkflowSubmitted)
      // Eventually the first workflow and its subworkflows should finish so the second workflow can run.
      waitForState(secondWorkflowId, WorkflowSucceeded)

      system.stop(rootActor)
    }
  }
}
