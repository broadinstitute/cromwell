package cromwell.backend.impl.jes

import java.util.UUID

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{ImplicitSender, TestDuration, TestKit}
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, ExecutionMode}
import cromwell.backend.async.ExecutionHandle
import cromwell.backend.impl.jes.JesAsyncBackendJobExecutionActor.JesPendingExecutionHandle
import cromwell.backend.impl.jes.RunStatus.Failed
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.prop.Tables.Table
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsValue}
import wdl4s.NamespaceWithWorkflow
import wdl4s.values.WdlString

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future, Promise}


class JesAsyncBackendJobExecutionActorSpec extends TestKit(ActorSystem("JesAsyncBackendJobExecutionActorSpec", ConfigFactory.parseString(
  """akka.loggers = ["akka.testkit.TestEventListener"]"""))) with WordSpecLike with Matchers with BeforeAndAfterAll with ImplicitSender with Mockito {

  val Timeout = 5.seconds.dilated

  val YoSup =
    """
      |task sup {
      |  String addressee
      |  command {
      |    echo "yo sup ${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |    [PREEMPTIBLE]
      |  }
      |}
      |
      |workflow sup {
      |  call sup
      |}
    """.stripMargin

  val Inputs = Map("sup.sup.addressee" -> WdlString("dog"))

  val NoOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

  val BackendConfig =
    """
      project = "my-cromwell-workflows"
      root = "gs://my-cromwell-workflows-bucket"

      genomics {
        auth = "application-default"
        endpoint-url = "https://genomics.googleapis.com/"
      }

      filesystems {
        gcs {
          auth = "application-default"
        }
      }
    """

  val AllConfig =
    """
       google {
         application-name = "cromwell"
         auths = [
           {
             name = "application-default"
             scheme = "application_default"
           }
         ]
       }

       backend {
         default = "JES"
         providers {
           JES {
             actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleFactory"
             config {""" + BackendConfig +
      """
             }
           }
         }
       }
      """.stripMargin

  val defaultBackendConfig = new BackendConfigurationDescriptor(ConfigFactory.parseString(BackendConfig), ConfigFactory.parseString(AllConfig))

  private def buildJobDescriptor(attempt: Int, preemptible: Int) = {
    val workflowDescriptor = new BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(YoSup.replace("[PREEMPTIBLE]", "preemptible: " + preemptible)),
      Inputs,
      NoOptions
    )

    val key = BackendJobDescriptorKey(workflowDescriptor.workflowNamespace.workflow.calls.head, None, attempt)
    BackendJobDescriptor(workflowDescriptor, key, Inputs)
  }

  private def executionActor(jobDescriptor: BackendJobDescriptor,
                             configurationDescriptor: BackendConfigurationDescriptor,
                             promise: Promise[BackendJobExecutionResponse],
                             errorCode: Int,
                             innerErrorCode: Int): ActorRef = {

    // Mock/stub out the bits that would reach out to JES.
    val run = mock[Run]
    run.status() returns new Failed(errorCode, Option(s"$innerErrorCode: I seen some things man"),
      Seq.empty)

    val handle = JesPendingExecutionHandle(jobDescriptor, Seq.empty, run, None)
    val jesConfiguration = new JesConfiguration(configurationDescriptor)
    class TestableJesJobExecutionActor extends JesAsyncBackendJobExecutionActor(jobDescriptor, promise, jesConfiguration) {
      // TODO: PBE: services are currently implemented in the engine, so we can't spin them up in specs
      override val serviceRegistryActor = system.actorOf(Props.empty)
      override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = Future.successful(handle)
    }
    system.actorOf(Props(new TestableJesJobExecutionActor), "TestableJesJobExecutionActor-" + UUID.randomUUID)
  }

  private def run(attempt: Int, preemptible: Int, errorCode: Int, innerErrorCode: Int): BackendJobExecutionResponse = {
    within(Timeout) {
      val promise = Promise[BackendJobExecutionResponse]()
      val jobDescriptor = buildJobDescriptor(attempt, preemptible)
      val backend = executionActor(jobDescriptor, defaultBackendConfig, promise, errorCode, innerErrorCode)
      backend ! Execute
      Await.result(promise.future, Duration.Inf)
    }
  }

  val expectations = Table(
    ("attempt", "preemptible", "errorCode", "innerErrorCode", "shouldRetry"),
    // No preemptible attempts allowed, nothing should be retryable.
    (       1,             0,          10,               13,        false),
    (       1,             0,          10,               14,        false),
    (       1,             0,          10,               15,        false),
    (       1,             0,          11,               13,        false),
    (       1,             0,          11,               14,        false),
    // 1 preemptible attempt allowed, but not all failures represent preemptions.
    (       1,             1,          10,               13,        true),
    (       1,             1,          10,               14,        true),
    (       1,             1,          10,               15,        false),
    (       1,             1,          11,               13,        false),
    (       1,             1,          11,               14,        false),
    // 1 preemptible attempt allowed, but now on the second attempt nothing should be retryable.
    (       2,             1,          10,               13,        false),
    (       2,             1,          10,               14,        false),
    (       2,             1,          10,               15,        false),
    (       2,             1,          11,               13,        false),
    (       2,             1,          11,               14,        false)
  )

  "JesAsyncBackendJobExecutionActor" should {
    "handle call failures appropriately with respect to preemption" in {
      expectations foreach { case (attempt, preemptible, errorCode, innerErrorCode, shouldRetry) =>
        run(attempt, preemptible, errorCode, innerErrorCode).getClass.getSimpleName match {
          case "FailedNonRetryableResponse" => shouldRetry shouldBe false
          case "FailedRetryableResponse" => shouldRetry shouldBe true
          case huh => fail(s"Unexpected response class name: '$huh'")
        }
      }
    }
  }
}
