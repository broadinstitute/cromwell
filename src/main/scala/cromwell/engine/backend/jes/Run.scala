package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.{CancelOperationRequest, Logging, RunPipelineRequest, ServiceAccount, _}
import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.jes.JesBackend.JesParameter
import cromwell.engine.backend.jes.Run.{Failed, Running, Success, _}
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.{JesCallBackendInfo, JesId, JesStatus}
import cromwell.engine.workflow.CallKey
import cromwell.util.google.GoogleScopes
import org.slf4j.LoggerFactory

import scala.annotation.tailrec
import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.language.postfixOps


object Run  {
  val JesServiceAccount = new ServiceAccount().setEmail("default").setScopes(GoogleScopes.Scopes.asJava)
  lazy val Log = LoggerFactory.getLogger("main")
  lazy val MaximumPollingInterval = ConfigFactory.load.getConfig("backend").getConfig("jes").getInt("maximumPollingInterval") * 1000
  val InitialPollingInterval = 5 seconds
  val PollingBackoffFactor = 1.1

  def apply(pipeline: Pipeline): Run = {
    val rpr = new RunPipelineRequest().setPipelineId(pipeline.id).setProjectId(pipeline.projectId).setServiceAccount(JesServiceAccount)
    val tag = s"JES Run [UUID(${pipeline.workflow.shortId}):${pipeline.key.scope.name}]"

    rpr.setInputs(pipeline.jesParameters.filter(_.isInput).toRunMap)
    Log.info(s"$tag Inputs:\n${stringifyMap(rpr.getInputs.asScala.toMap)}")

    rpr.setOutputs(pipeline.jesParameters.filter(_.isOutput).toRunMap)
    Log.info(s"$tag Outputs:\n${stringifyMap(rpr.getOutputs.asScala.toMap)}")

    val logging = new Logging()
    logging.setGcsPath(pipeline.gcsPath)
    rpr.setLogging(logging)

    // Currently, some resources (specifically disk) need to be specified both at pipeline creation and pipeline run time
    rpr.setResources(pipeline.runtimeInfo.resources)

    val id = pipeline.genomicsService.pipelines().run(rpr).execute().getName
    Log.info(s"$tag JES ID is $id")
    new Run(id, pipeline, tag)
  }

  private def stringifyMap(m: Map[String, String]): String = m map { case(k, v) => s"  $k -> $v"} mkString "\n"

  implicit class RunJesParameters(val params: Seq[JesParameter]) extends AnyVal {
    def toRunMap = params.map(p => p.name -> p.gcs).toMap.asJava
  }

  implicit class RunOperationExtension(val operation: Operation) extends AnyVal {
    def hasStarted = operation.getMetadata.asScala.get("startTime") isDefined
  }

  sealed trait RunStatus
  trait TerminalRunStatus extends RunStatus
  case object Initializing extends RunStatus
  case object Running extends RunStatus
  case object Success extends TerminalRunStatus
  final case class Failed(errorCode: Int, errorMessage: String) extends TerminalRunStatus {
    // Don't want to include errorMessage or code in the snappy status toString:
    override def toString = "Failed"
  }

  @tailrec
  private final def waitForStatus(run: Run, previousStatus: Option[RunStatus], pollingInterval: Double, breakout: RunStatus => Boolean): RunStatus = {
    val currentStatus = run.status()

    if (!(previousStatus contains currentStatus)) {
      // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise
      // just use the state names.
      val prevStateName = previousStatus map { _.toString } getOrElse "-"
      Log.info(s"${run.tag}: Status change from $prevStateName to $currentStatus")

      // Update the database state:
      val newBackendInfo = JesCallBackendInfo(Option(JesId(run.jesId)), Option(JesStatus(currentStatus.toString)))
      globalDataAccess.updateExecutionBackendInfo(run.workflowId, CallKey(run.call, run.pipeline.key.index), newBackendInfo)
    }

    if (breakout(currentStatus)) {
      currentStatus
    } else {
      Thread.sleep(pollingInterval.toInt)
      waitForStatus(run, Option(currentStatus), nextPollingInterval(pollingInterval, MaximumPollingInterval), breakout)
    }
  }

  private final def nextPollingInterval(previousPollingInterval: Double, maximumPollingInterval: Int): Double = {
    Math.min(previousPollingInterval * PollingBackoffFactor, maximumPollingInterval)
  }
}

case class Run(jesId: String, pipeline: Pipeline, tag: String) {

  lazy val workflowId = pipeline.workflow.id
  lazy val call = pipeline.key.scope

  def status(): RunStatus = {
    val op = pipeline.genomicsService.operations().get(jesId).execute

    if (op.getDone) {
      // If there's an error, generate a Failed status. Otherwise, we were successful!
      Option(op.getError) map { x => Failed(x.getCode, x.getMessage) } getOrElse Success
    } else if (op.hasStarted) {
      Running
    } else {
      Initializing
    }
  }

  final def waitUntilComplete(previousStatus: RunStatus): TerminalRunStatus = {
    val terminalStatus = Run.waitForStatus(this, Option(previousStatus), InitialPollingInterval.toMillis, {
      case x: TerminalRunStatus => true
      case _ => false
    })

    terminalStatus match {
      case x: TerminalRunStatus => x
      case _ => Failed(-1, "Unexpectedly stopped checking status.") // Assuming waitForStatus works, this never happens.
    }
  }

  final def waitUntilRunningOrComplete: RunStatus = Run.waitForStatus(this, None, InitialPollingInterval.toMillis, {
    case Running => true
    case x: TerminalRunStatus => true
    case _ => false
  })

  def abort(): Unit = {
    val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
    pipeline.genomicsService.operations().cancel(jesId, cancellationRequest).execute
  }
}