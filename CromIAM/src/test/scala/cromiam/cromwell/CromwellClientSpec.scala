package cromiam.cromwell

import java.net.URL

import akka.actor.{ActorRef, ActorSystem}
import akka.event.NoLogging
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken}
import akka.http.scaladsl.model.{HttpHeader, HttpRequest, HttpResponse}
import akka.stream.ActorMaterializer
import cats.data.NonEmptyList
import cats.effect.IO
import cromiam.auth.User
import cromwell.api.model._
import cromwell.api.{CromwellClient => CromwellApiClient}
import cromwell.services.instrumentation.CromwellInstrumentation.InstrumentationPath
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId
import org.scalatest.{AsyncFlatSpec, BeforeAndAfterAll, Matchers}
import spray.json.{JsObject, JsString}

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, ExecutionContextExecutor}

class CromwellClientSpec extends AsyncFlatSpec with Matchers with BeforeAndAfterAll {
  import CromwellClientSpec._

  implicit val actorSystem: ActorSystem = ActorSystem("CromwellClientSpec")
  implicit val ece: ExecutionContextExecutor = actorSystem.dispatcher
  implicit val materializer: ActorMaterializer = ActorMaterializer()

  val cromwellClient = new MockCromwellClient()

  val fakeHttpRequest = HttpRequest()

  override protected def afterAll(): Unit = {
    actorSystem.terminate()
    super.afterAll()
  }

  "CromwellClient" should "eventually return a subworkflow's root workflow id" in {
    cromwellClient.getRootWorkflow(SubworkflowId.id.toString, FictitiousUser, fakeHttpRequest).map(w => assert(w == RootWorkflowId.id.toString))
      .asIo.unsafeToFuture()
  }

  it should "eventually return a top level workflow's ID when requesting root workflow id" in {
    cromwellClient.getRootWorkflow(RootWorkflowId.id.toString, FictitiousUser, fakeHttpRequest).map(w => assert(w == RootWorkflowId.id.toString))
      .asIo.unsafeToFuture()
  }

  it should "properly fetch the collection for a workflow with a collection name" in {
    cromwellClient.collectionForWorkflow(RootWorkflowId.id.toString, FictitiousUser, fakeHttpRequest).map(c =>
      assert(c.name == CollectionName)
    ).asIo.unsafeToFuture()
  }

  it should "throw an exception if the workflow doesn't have a collection" in {
    recoverToExceptionIf[IllegalArgumentException] {
      cromwellClient.collectionForWorkflow(WorkflowIdWithoutCollection.id.toString, FictitiousUser, fakeHttpRequest)
        .asIo.unsafeToFuture()
    } map { exception =>
      assert(exception.getMessage == s"Workflow $WorkflowIdWithoutCollection has no associated collection")
    }
  }
}

object CromwellClientSpec {
  final class MockCromwellClient()(implicit system: ActorSystem,
                                   ece: ExecutionContextExecutor,
                                   materializer: ActorMaterializer)
  extends CromwellClient("http", "bar", 1, NoLogging, ActorRef.noSender) {
    override val cromwellApiClient: CromwellApiClient = new MockCromwellApiClient()

    override def sendTimingApi(statsDPath: InstrumentationPath,
                               timing: FiniteDuration,
                               prefixToStatsd: NonEmptyList[String]
                              ): Unit = ()
  }

  final class MockCromwellApiClient()(implicit actorSystem: ActorSystem, materializer: ActorMaterializer)
    extends CromwellApiClient(new URL("http://foo.com"), "bar") {


    override def labels(workflowId: WorkflowId, headers: List[HttpHeader] = defaultHeaders)
                       (implicit ec: ExecutionContext): FailureResponseOrT[WorkflowLabels] = {
      if (workflowId == RootWorkflowId) {
        FailureResponseOrT.pure(FictitiousWorkflowLabelsWithCollection)
      } else if (workflowId == WorkflowIdWithoutCollection) {
        FailureResponseOrT.pure(FictitiousWorkflowLabelsWithoutCollection)
      } else {
        FailureResponseOrT[IO, HttpResponse, WorkflowLabels] {
          IO.raiseError(new RuntimeException("Unexpected workflow ID sent to MockCromwellApiClient"))
        }
      }
    }

    override def metadata(workflowId: WorkflowId,
                 args: Option[Map[String, List[String]]] = None,
                 headers: List[HttpHeader] = defaultHeaders
                )(implicit ec: ExecutionContext): FailureResponseOrT[WorkflowMetadata] = {
      if (workflowId == RootWorkflowId) FailureResponseOrT.pure(RootWorkflowMetadata)
      else if (workflowId == SubworkflowId) FailureResponseOrT.pure(SubWorkflowMetadata)
      else FailureResponseOrT[IO, HttpResponse, WorkflowMetadata] {
        IO.raiseError(new RuntimeException("Unexpected workflow ID sent to MockCromwellApiClient"))
      }
    }
  }

  private val SubworkflowId = WorkflowId.fromString("58114f5c-f439-4488-8d73-092273cf92d9")
  private val RootWorkflowId = WorkflowId.fromString("998d137e-7213-44ac-8f6f-24e6e23adaa5")
  private val WorkflowIdWithoutCollection = WorkflowId.fromString("472234d-7213-44ac-8f6f-24e6e23adaa5")

  val FictitiousUser = User(WorkbenchUserId("123456780"), Authorization(OAuth2BearerToken("my-token")))

  val RootWorkflowMetadata = WorkflowMetadata("""{
                                               "workflowName": "wf_echo",
                                               "calls": {},
                                               "id": "998d137e-7213-44ac-8f6f-24e6e23adaa5"
                                             }""")

  val SubWorkflowMetadata = WorkflowMetadata("""{
                                                "workflowName": "wf_echo",
                                                "calls": {},
                                                "id": "58114f5c-f439-4488-8d73-092273cf92d9",
                                                "rootWorkflowId": "998d137e-7213-44ac-8f6f-24e6e23adaa5"
                                              }""")

  val CollectionName = "foo"
  val FictitiousWorkflowLabelsWithCollection = WorkflowLabels(RootWorkflowId.id.toString, JsObject(Map("caas-collection-name" -> JsString(CollectionName))))
  val FictitiousWorkflowLabelsWithoutCollection = WorkflowLabels(RootWorkflowId.id.toString, JsObject(Map("something" -> JsString("foo"))))
}

