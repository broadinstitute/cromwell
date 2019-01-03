package cromiam.cromwell

import java.net.URL

import akka.actor.ActorSystem
import akka.event.NoLogging
import akka.http.scaladsl.model.HttpHeader
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken}
import akka.stream.ActorMaterializer
import cromiam.auth.User
import cromwell.api.model.{WorkflowId, WorkflowLabels, WorkflowMetadata}
import org.scalatest.{AsyncFlatSpec, BeforeAndAfterAll, Matchers}
import cromwell.api.{CromwellClient => CromwellApiClient}
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId
import spray.json.{JsObject, JsString}

import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future}

class CromwellClientSpec extends AsyncFlatSpec with Matchers with BeforeAndAfterAll {
  import CromwellClientSpec._

  import ExecutionContext.Implicits.global

  implicit val actorSystem: ActorSystem = ActorSystem("CromwellClientSpec")
  implicit val ece: ExecutionContextExecutor = actorSystem.dispatcher
  implicit val materializer: ActorMaterializer = ActorMaterializer()

  val cromwellClient = new MockCromwellClient()

  override protected def afterAll(): Unit = {
    actorSystem.terminate()
    super.afterAll()
  }

  "CromwellClient" should "eventually return a subworkflow's root workflow id" in {
    cromwellClient.getRootWorkflow(SubworkflowId.id.toString, FictitiousUser).map(w => assert(w == RootWorkflowId.id.toString))
  }

  it should "eventually return a top level workflow's ID when requesting root workflow id" in {
    cromwellClient.getRootWorkflow(RootWorkflowId.id.toString, FictitiousUser).map(w => assert(w == RootWorkflowId.id.toString))
  }

  it should "properly fetch the collection for a workflow with a collection name" in {
    cromwellClient.collectionForWorkflow(RootWorkflowId.id.toString, FictitiousUser).map(c =>
      assert(c.name == CollectionName)
    )
  }

  it should "throw an exception if the workflow doesn't have a collection" in {
    recoverToExceptionIf[IllegalArgumentException] {
      cromwellClient.collectionForWorkflow(WorkflowIdWithoutCollection.id.toString, FictitiousUser)
    } map { exception =>
      assert(exception.getMessage == s"Workflow $WorkflowIdWithoutCollection has no associated collection")
    }
  }
}

object CromwellClientSpec {
  final class MockCromwellClient()(implicit system: ActorSystem,
                                   ece: ExecutionContextExecutor,
                                   materializer: ActorMaterializer)
  extends CromwellClient("http", "bar", 1, NoLogging) {
    override val cromwellApiClient: CromwellApiClient = new MockCromwellApiClient()
  }

  final class MockCromwellApiClient()(implicit actorSystem: ActorSystem, materializer: ActorMaterializer)
    extends CromwellApiClient(new URL("http://foo.com"), "bar") {

    override def labels(workflowId: WorkflowId, headers: List[HttpHeader] = defaultHeaders)(implicit ec: ExecutionContext): Future[WorkflowLabels] = {
      if (workflowId == RootWorkflowId) Future.successful(FictitiousWorkflowLabelsWithCollection)
      else if (workflowId == WorkflowIdWithoutCollection) Future.successful(FictitiousWorkflowLabelsWithoutCollection)
      else Future.failed(new RuntimeException("Unexpected workflow ID sent to MockCromwellApiClient"))
    }

    override def metadata(workflowId: WorkflowId,
                 args: Option[Map[String, List[String]]] = None,
                 headers: List[HttpHeader] = defaultHeaders
                )(implicit ec: ExecutionContext): Future[WorkflowMetadata] = {
      if (workflowId == RootWorkflowId) Future.successful(RootWorkflowMetadata)
      else if (workflowId == SubworkflowId) Future.successful(SubWorkflowMetadata)
      else Future.failed(new RuntimeException("Unexpected workflow ID sent to MockCromwellApiClient"))
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

