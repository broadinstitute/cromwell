package cromwell.backend.impl.spark

import java.nio.file.Path

import akka.testkit.ImplicitSender
import better.files._
import Cmds._

import scala.concurrent.{Future, Promise}
import cromwell.core.TestKitSuite
import org.mockito.Matchers._
import org.mockito.Mockito
import org.mockito.Mockito._
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.mockito.MockitoSugar
import org.scalatest.{BeforeAndAfter, Matchers, WordSpecLike}
import spray.json._

import scala.concurrent.duration._
import cromwell.backend.impl.spark.SparkClusterProcess.{Failed, _}
import org.scalatest.concurrent.ScalaFutures
import spray.http._
import SparkClusterJsonProtocol._

class SparkClusterProcessSpec extends TestKitSuite("SparkClusterProcess")
  with WordSpecLike
  with Matchers
  with MockitoSugar
  with BeforeAndAfter
  with ImplicitSender
  with ScalaFutures {

  private val timeout = Timeout(5.seconds)

  var testRoot: File = _
  var jsonFile: File = _
  var rcFile: File = _

  val jsonFileName = "test.json"
  val rcFileName = "rc"

  before {
    testRoot = File.newTemporaryDirectory("spark-test")
    jsonFile = testRoot / jsonFileName
    rcFile = testRoot / rcFileName
  }

  after {
    rm(testRoot)
  }

  private val sampleSubmissionResponse =
    """
      |Running Spark using the REST application submission protocol.
      |16/08/06 18:35:26 INFO rest.RestSubmissionClient: Submitting a request to launch an application in spark://host-10-0-1-53:6066.
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: Submission successfully created as driver-20160806183527-0006. Polling submission state...
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: Submitting a request for the status of submission driver-20160806183527-0006 in spark://host-10-0-1-53:6066.
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: State of driver driver-20160806183527-0006 is now RUNNING.
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: Driver is running on worker worker-20160801162431-10.0.1.55-43834 at 10.0.1.55:43834.
      |16/08/06 18:35:27 INFO rest.RestSubmissionClient: Server responded with CreateSubmissionResponse:
      |{
      |  "action" : "CreateSubmissionResponse",
      |  "message" : "Driver successfully submitted as driver-20160806183527-0006",
      |  "serverSparkVersion" : "1.6.1",
      |  "submissionId" : "driver-20160806183527-0006",
      |  "success" : true
      |}
    """.stripMargin

  private val mockSuccessClusterResponse = SparkDriverStateQueryResponse(action = "SubmissionStatusResponse", driverState = "FINISHED", serverSparkVersion = "1.6.1",
    submissionId = "driver-20160803181054-0000", success = true, workerHostPort = "10.0.1.55:43834", workerId = "worker-20160801162431-10.0.1.55-43834")
  private val mockFailedClusterResponse = SparkDriverStateQueryResponse(action = "SubmissionStatusResponse", driverState = "FAILED", serverSparkVersion = "1.6.1",
    submissionId = "driver-20160803181054-0000", success = true, workerHostPort = "10.0.1.55:43834", workerId = "worker-20160801162431-10.0.1.55-43834")
  private val mockUnknownResponse = SparkDriverStateQueryResponse(action = "SubmissionStatusResponse", driverState = "Unknown", serverSparkVersion = "1.6.1",
    submissionId = "driver-20160803181054-0000", success = true, workerHostPort = "10.0.1.55:43834", workerId = "worker-20160801162431-10.0.1.55-43834")
  private val mockRunningClusterResponse = SparkDriverStateQueryResponse(action = "SubmissionStatusResponse", driverState = "RUNNING", serverSparkVersion = "1.6.1",
    submissionId = "driver-20160803181054-0000", success = true, workerHostPort = "10.0.1.55:43834", workerId = "worker-20160801162431-10.0.1.55-43834")

  private val mockSuccessHttpResponse = HttpResponse(StatusCodes.OK, HttpEntity(ContentTypes.`application/json`, mockSuccessClusterResponse.toJson.toString))
  private val mockRunningHttpResponse = HttpResponse(StatusCodes.OK, HttpEntity(ContentTypes.`application/json`, mockRunningClusterResponse.toJson.toString))
  private val mockFailedHttpResponse = HttpResponse(StatusCodes.OK, HttpEntity(ContentTypes.`application/json`, mockFailedClusterResponse.toJson.toString))
  private val mockBadHttpResponse = HttpResponse(StatusCodes.BadRequest)

  "SparkCluster process" should {
    "return finished status when Spray rest client return successful http response" in {
      val sparkClusterProcess = Mockito.spy(new SparkClusterProcess()(system))
      jsonFile write sampleSubmissionResponse
      doReturn(Future.successful(mockSuccessHttpResponse)).when(sparkClusterProcess).makeHttpRequest(any[HttpRequest])
      whenReady(sparkClusterProcess.startMonitoringSparkClusterJob(testRoot.path, jsonFileName), timeout) { response =>
        response shouldBe Finished
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(testRoot.path, jsonFileName)
        verify(sparkClusterProcess, times(1)).parseJsonForSubmissionIdAndStatus(any[Path])
        verify(sparkClusterProcess, times(1)).monitorSparkClusterJob(any[String], any[Path], any[Promise[Unit]])
        verify(sparkClusterProcess, times(1)).pollForJobStatus(any[String])
      }
    }

    "return failed status when Spray rest client return unsuccessful http response" in {
      val sparkClusterProcess = Mockito.spy(new SparkClusterProcess()(system))
      jsonFile write sampleSubmissionResponse
      doReturn(Future.successful(mockFailedHttpResponse)).when(sparkClusterProcess).makeHttpRequest(any[HttpRequest])
      whenReady(sparkClusterProcess.startMonitoringSparkClusterJob(testRoot.path, jsonFileName), timeout) { response =>
        response shouldBe a[Failed]
        assert(response.asInstanceOf[Failed].error.getMessage.contains("Spark Driver returned failed status"))
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(testRoot.path, jsonFileName)
        verify(sparkClusterProcess, times(1)).parseJsonForSubmissionIdAndStatus(any[Path])
        verify(sparkClusterProcess, times(1)).monitorSparkClusterJob(any[String], any[Path], any[Promise[Unit]])
        verify(sparkClusterProcess, times(1)).pollForJobStatus(any[String])
      }
    }

    "return failed status when Spray rest client return bad http response" in {
      val sparkClusterProcess = Mockito.spy(new SparkClusterProcess()(system))
      jsonFile write sampleSubmissionResponse
      doReturn(Future.successful(mockBadHttpResponse)).when(sparkClusterProcess).makeHttpRequest(any[HttpRequest])
      whenReady(sparkClusterProcess.startMonitoringSparkClusterJob(testRoot.path, jsonFileName), timeout) { response =>
        response shouldBe a[Failed]
        assert(response.asInstanceOf[Failed].error.getMessage.contains("Spark Driver returned failed status"))
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(testRoot.path, jsonFileName)
        verify(sparkClusterProcess, times(1)).parseJsonForSubmissionIdAndStatus(any[Path])
        verify(sparkClusterProcess, times(1)).monitorSparkClusterJob(any[String], any[Path], any[Promise[Unit]])
        verify(sparkClusterProcess, times(1)).pollForJobStatus(any[String])
      }
    }

    "return failed status when makeHttpRequest throws exception" in {
      val sparkClusterProcess = Mockito.spy(new SparkClusterProcess()(system))
      jsonFile write sampleSubmissionResponse
      doReturn(Future.failed(new IllegalStateException("request timed out"))).when(sparkClusterProcess).makeHttpRequest(any[HttpRequest])
      whenReady(sparkClusterProcess.startMonitoringSparkClusterJob(testRoot.path, jsonFileName), timeout) { response =>
        response shouldBe a[Failed]
        assert(response.asInstanceOf[Failed].error.getMessage.contains("Spark Driver returned failed status"))
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(testRoot.path, jsonFileName)
        verify(sparkClusterProcess, times(1)).parseJsonForSubmissionIdAndStatus(any[Path])
        verify(sparkClusterProcess, times(1)).monitorSparkClusterJob(any[String], any[Path], any[Promise[Unit]])
        verify(sparkClusterProcess, times(1)).pollForJobStatus(any[String])
      }
    }

    "return failed status when parser failed" in {
      val sparkClusterProcess = Mockito.spy(new SparkClusterProcess()(system))
      jsonFile write ""
      whenReady(sparkClusterProcess.startMonitoringSparkClusterJob(testRoot.path, jsonFileName), timeout) { response =>
        response shouldBe a[Failed]
        assert(response.asInstanceOf[Failed].error.getMessage.contains("Unable to get json out of submission response file"))
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(testRoot.path, jsonFileName)
        verify(sparkClusterProcess, times(1)).parseJsonForSubmissionIdAndStatus(any[Path])
      }
    }

    "return failed status when Spray rest client return unknown response" in {
      val sparkClusterProcess = Mockito.spy(new SparkClusterProcess()(system))
      jsonFile write sampleSubmissionResponse
      doReturn(Future.successful(mockUnknownResponse)).when(sparkClusterProcess).makeHttpRequest(any[HttpRequest])
      whenReady(sparkClusterProcess.startMonitoringSparkClusterJob(testRoot.path, jsonFileName), timeout) { response =>
        response shouldBe a[Failed]
        assert(response.asInstanceOf[Failed].error.getMessage.contains("Spark Driver returned failed status"))
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(testRoot.path, jsonFileName)
        verify(sparkClusterProcess, times(1)).parseJsonForSubmissionIdAndStatus(any[Path])
        verify(sparkClusterProcess, times(1)).monitorSparkClusterJob(any[String], any[Path], any[Promise[Unit]])
        verify(sparkClusterProcess, times(1)).pollForJobStatus(any[String])
      }
    }

    "return finished status when Spray rest client return running and then finished response" in {
      val sparkClusterProcess = Mockito.spy(new SparkClusterProcess()(system))
      jsonFile write sampleSubmissionResponse
      doReturn(Future.successful(mockRunningHttpResponse)).doReturn(Future.successful(mockSuccessHttpResponse))
        .when(sparkClusterProcess).makeHttpRequest(any[HttpRequest])
      whenReady(sparkClusterProcess.startMonitoringSparkClusterJob(testRoot.path, jsonFileName), timeout) { response =>
        response shouldBe Finished
        verify(sparkClusterProcess, times(1)).startMonitoringSparkClusterJob(testRoot.path, jsonFileName)
        verify(sparkClusterProcess, times(1)).parseJsonForSubmissionIdAndStatus(any[Path])
        verify(sparkClusterProcess, times(2)).monitorSparkClusterJob(any[String], any[Path], any[Promise[Unit]])
        verify(sparkClusterProcess, times(2)).pollForJobStatus(any[String])
      }
    }

  }

}
