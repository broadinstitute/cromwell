package cromwell.backend.google.pipelines.v1alpha2.api

import java.net.URL

import akka.actor.ActorRef
import akka.testkit.{TestActorRef, TestProbe}
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.services.genomics.model.Operation
import cromwell.backend.google.pipelines.common.PipelinesApiTestConfig.papiConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIStatusPollRequest
import cromwell.backend.google.pipelines.common.api.{PipelinesApiRequestManager, PipelinesApiRequestWorkerSpec, TestPipelinesApiBatchHandler, TestPipelinesApiRequestWorker}
import cromwell.backend.google.pipelines.v1alpha2.api.request.RequestHandler

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Success, Try}

class PipelinesApiRequestWorkerV1Spec extends PipelinesApiRequestWorkerSpec[Operation] {
  override implicit var batchHandler: TestPipelinesApiBatchHandler[Operation] = _

  before {
    managerProbe = TestProbe()
    val registryProbe = TestProbe()
    batchHandler = new TestPipelinesApiBatchHandlerV1(managerProbe.ref)
    workerActor = TestActorRef(TestPipelinesApiRequestWorker.props(managerProbe.ref, papiConfiguration, registryProbe.ref)(batchHandler), managerProbe.ref)
  }
}

class TestPipelinesApiBatchHandlerV1(pollingManager: ActorRef) extends TestPipelinesApiBatchHandler[Operation] {
  val delegate = new RequestHandler("cromwell-test-app", new URL("https://genomics.googleapis.com/")) {
    override def interpretOperationStatus(operation: Operation) = {
      mockStatusInterpreter(operation)
    }
    override protected def mkErrorString(e: GoogleJsonError) = "N/A"
  }
  
  override def statusPollResultHandler(pollRequest: PipelinesApiRequestManager.PAPIStatusPollRequest, completionPromise: Promise[Try[Unit]]) = {
    delegate.statusPollResultHandler(pollRequest, completionPromise, pollingManager)
  }
  
  override def makeBatchRequest = null
  override def enqueue[T <: PipelinesApiRequestManager.PAPIApiRequest](papiApiRequest: T, batchRequest: BatchRequest, pollingManager: ActorRef)
                                                                      (implicit ec: ExecutionContext)= {
    papiApiRequest match {
      case poll: PAPIStatusPollRequest => enqueueStatusPollInBatch(poll, batchRequest)
      case _ => Future.successful(Success(()))
    }
  }
}
