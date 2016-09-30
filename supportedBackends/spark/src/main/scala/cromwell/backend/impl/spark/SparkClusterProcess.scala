package cromwell.backend.impl.spark

import java.nio.file.Path

import akka.actor.ActorSystem
import cromwell.backend.impl.spark.SparkClusterProcess.{SparkJobSubmissionResponse, TerminalStatus}
import spray.http.{HttpRequest, HttpResponse, StatusCodes}
import spray.json.{DefaultJsonProtocol, JsonParser}
import spray.client.pipelining._

import scala.concurrent.{ExecutionContext, Future, Promise}
import better.files._
import com.typesafe.scalalogging.Logger
import org.slf4j.LoggerFactory
import scala.concurrent.duration._
import scala.util.{Failure, Success, Try}

object SparkClusterProcess {

  case class SparkDriverStateQueryResponse(action: String, driverState: String, serverSparkVersion: String,
                                           submissionId: String, success: Boolean, workerHostPort: String, workerId: String)

  case class SparkJobSubmissionResponse(action: String, message: String, serverSparkVersion: String, submissionId: String, success: Boolean)

  sealed trait TerminalStatus

  case class Failed(error: Throwable) extends TerminalStatus

  case object Finished extends TerminalStatus

  object SparkClusterJsonProtocol extends DefaultJsonProtocol {
    implicit val sparkStatusResponseFormat = jsonFormat7(SparkDriverStateQueryResponse)
    implicit val successfulParsedResponseFormat = jsonFormat5(SparkJobSubmissionResponse)
  }

}

trait SparkClusterRestClient {
  def sendAndReceive: SendReceive

  def makeHttpRequest(httpRequest: HttpRequest): Future[HttpResponse]
}

trait SparkClusterProcessMonitor {
  def startMonitoringSparkClusterJob(jobPath: Path, jsonFile: String): Future[TerminalStatus]

  def monitorSparkClusterJob(subId: String, rcPath: Path, promise: Promise[Unit]): Unit

  def evaluateMonitoringFuture(rcPath: Path): Unit

  def completeMonitoringProcess(rcPath: Path, status: String, promise: Promise[Unit]): Unit
}

trait SparkClusterJobParser {
  def parseJsonForSubmissionIdAndStatus(jsonFile: Path): SparkJobSubmissionResponse
}

class SparkClusterProcess(implicit system: ActorSystem) extends SparkProcess
  with SparkClusterRestClient with SparkClusterJobParser with SparkClusterProcessMonitor {

  import SparkClusterProcess._
  import spray.httpx.SprayJsonSupport._
  import SparkClusterJsonProtocol._

  implicit lazy val ec: ExecutionContext = system.dispatcher
  lazy val completionPromise = Promise[TerminalStatus]()
  lazy val monitorPromise = Promise[Unit]()
  val tag = this.getClass.getSimpleName
  lazy val logger = Logger(LoggerFactory.getLogger(getClass.getName))

  override def sendAndReceive: SendReceive = sendReceive

  override def startMonitoringSparkClusterJob(jobPath: Path, jsonFile: String): Future[TerminalStatus] = {
    Future(parseJsonForSubmissionIdAndStatus(jobPath.resolve(jsonFile))) onComplete {
      case Success(resp: SparkJobSubmissionResponse) =>
        val rcPath = jobPath.resolve("rc")
        monitorSparkClusterJob(resp.submissionId, rcPath, monitorPromise)
        evaluateMonitoringFuture(rcPath)
      case Failure(exception: Throwable) =>
        logger.error(s"{} Spark Job failed to submit successfully following reason: {}", tag, exception.getMessage)
        completionPromise success Failed(exception)
    }
    completionPromise.future
  }

  override def monitorSparkClusterJob(subId: String, rcPath: Path, promise: Promise[Unit]): Unit = {
    pollForJobStatus(subId) onComplete {
      case Success(resp: SparkDriverStateQueryResponse) if resp.driverState == "RUNNING" =>
        logger.debug(s"{} Spark Driver is now in :{} state", tag, "RUNNING")
        system.scheduler.scheduleOnce(1000.milliseconds) {
          monitorSparkClusterJob(subId, rcPath, promise)
        }
      case Success(resp: SparkDriverStateQueryResponse) if resp.driverState == "FINISHED" =>
        logger.debug(s"{} Spark Driver is now in :{} state", tag, "FINISHED")
        completeMonitoringProcess(rcPath, "0", promise)
      case Success(resp: SparkDriverStateQueryResponse) if resp.driverState == "FAILED" =>
        logger.debug(s"{} Spark Driver is now in :{} state", tag, "FAILED")
        completeMonitoringProcess(rcPath, "-1", promise)
      case Success(_) =>
        logger.error(s"{} Spark Driver is now in Unknown state as polling returned success", tag)
        completeMonitoringProcess(rcPath, "-1", promise)
      case Failure(error) =>
        logger.error(s"{} Spark cluster poll for status failed Reason: {} and error: {}", tag, error.getMessage, error)
        completeMonitoringProcess(rcPath, "-1", promise)
    }
  }

  override def evaluateMonitoringFuture(rcPath: Path) = {
    monitorPromise.future.onComplete {
      case Success(_) if File(rcPath).contentAsString.stripLineEnd.toInt == 0 =>
        completionPromise success Finished
      case Success(_) =>
        completionPromise success Failed(new IllegalStateException("Spark Driver returned failed status"))
      case Failure(err) =>
        completionPromise success Failed(err)
    }
  }

  override def completeMonitoringProcess(rcPath: Path, status: String, promise: Promise[Unit]) = {
    File(rcPath) write status
    val unitValue = ()
    promise success unitValue
    ()
  }

  def pollForJobStatus(subId: String): Future[SparkDriverStateQueryResponse] = {
    // on failure we use spark-master as a default value for the master hostname
    val sparkClusterMasterHostName = Try(sys.env("HOSTNAME")) match {
      case Success(s) => Option(s)
      case Failure(_) => None
    }

    val request = sparkClusterMasterHostName match {
      case Some(master) =>
        Get(s"http://$master:6066/v1/submissions/status/$subId")
      case None =>
        Get(s"http://spark-master:6066/v1/submissions/status/$subId")
    }

    makeHttpRequest(request) flatMap { v =>
      v.status match {
        case StatusCodes.OK => Future(v ~> unmarshal[SparkDriverStateQueryResponse])
        case _ =>
          val msg = s"Unexpected response received in response from Spark rest api. Response: $v"
          logger.error("{} reason: {}", tag, msg)
          throw new IllegalStateException(msg)
      }
    } recover {
      case error => throw new IllegalStateException(s"Reason: ${error.getMessage}", error)
    }
  }

  override def parseJsonForSubmissionIdAndStatus(jsonFile: Path): SparkJobSubmissionResponse = {
    val lines = File(jsonFile).contentAsString
    val sparkClusterSubmissionResponseRegex = """(?s)\{(.*)}""".r.unanchored
    val line = sparkClusterSubmissionResponseRegex findFirstIn lines match {
      case Some(content) => content
      case None =>
        val msg = "Unable to get json out of submission response file"
        logger.error("{} reason: {}", tag, msg)
        throw new IllegalStateException(msg)
    }

    JsonParser(line).convertTo[SparkJobSubmissionResponse]
  }

  override def makeHttpRequest(httpRequest: HttpRequest): Future[HttpResponse] = {
    val headers = httpRequest.headers
    sendAndReceive(httpRequest.withHeaders(headers))
  }
}
