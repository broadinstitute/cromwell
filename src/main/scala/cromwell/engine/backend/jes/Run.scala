package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.{CancelOperationRequest, Logging, RunPipelineRequest, ServiceAccount, _}
import com.typesafe.config.ConfigFactory
import cromwell.engine.AbortFunction
import cromwell.engine.backend.jes.JesBackend.{JesInput, JesOutput, JesParameter}
import cromwell.engine.backend.jes.Run.{Failed, Running, Success, _}
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.{JesCallBackendInfo, JesId, JesStatus}
import cromwell.engine.workflow.CallKey
import cromwell.util.google.GoogleScopes
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.language.postfixOps

object Run  {
  val JesServiceAccount = new ServiceAccount().setEmail("default").setScopes(GoogleScopes.Scopes.asJava)
  lazy val Log = LoggerFactory.getLogger("main")
  lazy val MaximumPollingInterval = Duration(ConfigFactory.load.getConfig("backend").getConfig("jes").getInt("maximumPollingInterval"), "seconds")
  val InitialPollingInterval = 5 seconds
  val PollingBackoffFactor = 1.1

  def apply(pipeline: Pipeline): Run = {
    val tag = s"JES Run [UUID(${pipeline.workflow.shortId}):${pipeline.key.tag}]"
    if (pipeline.pipelineId.isDefined == pipeline.runIdForResumption.isDefined) {
      val message =
        s"""
          |$tag: Exactly one of JES pipeline ID or run ID for resumption must be specified to create a Run.
          |pipelineId = ${pipeline.pipelineId}, runIdForResumption = ${pipeline.runIdForResumption}.
        """.stripMargin
      throw new RuntimeException(message)
    }

    def runPipeline: String = {

      val rpr = new RunPipelineRequest().setPipelineId(pipeline.pipelineId.get).setProjectId(pipeline.projectId).setServiceAccount(JesServiceAccount)

      rpr.setInputs(pipeline.jesParameters.collect({ case i: JesInput => i }).toRunMap)
      Log.info(s"$tag Inputs:\n${stringifyMap(rpr.getInputs.asScala.toMap)}")

      rpr.setOutputs(pipeline.jesParameters.collect({ case i: JesOutput => i }).toRunMap)
      Log.info(s"$tag Outputs:\n${stringifyMap(rpr.getOutputs.asScala.toMap)}")

      val logging = new Logging()
      logging.setGcsPath(pipeline.gcsPath)
      rpr.setLogging(logging)

      // Currently, some resources (specifically disk) need to be specified both at pipeline creation and pipeline run time
      rpr.setResources(pipeline.runtimeInfo.resources)

      val runId = pipeline.genomicsService.pipelines().run(rpr).execute().getName
      Log.info(s"$tag JES run ID is $runId")
      runId
    }

    // Only run the pipeline if the pipeline ID is defined.  The pipeline ID not being defined corresponds to a
    // resumption of a previous run, and runIdForResumption will be defined.  The Run code takes care of polling
    // in both the newly created and resumed scenarios.
    val runId = if (pipeline.pipelineId.isDefined) runPipeline else pipeline.runIdForResumption.get
    new Run(runId, pipeline, tag)
  }

  private def stringifyMap(m: Map[String, String]): String = m map { case(k, v) => s"  $k -> $v"} mkString "\n"

  implicit class RunJesParameters(val params: Seq[JesParameter]) extends AnyVal {
    def toRunMap = params.map(p => p.name -> p.gcs).toMap.asJava
  }

  implicit class RunOperationExtension(val operation: Operation) extends AnyVal {
    def hasStarted = operation.getMetadata.asScala.get("startTime") isDefined
  }

  sealed trait RunStatus {
    // Could be defined as false for Initializing and true otherwise, but this is more defensive.
    def isRunningOrComplete = this match {
      case Running | _: TerminalRunStatus => true
      case _ => false
    }
  }
  trait TerminalRunStatus extends RunStatus
  case object Initializing extends RunStatus
  case object Running extends RunStatus
  case object Success extends TerminalRunStatus
  final case class Failed(errorCode: Int, errorMessage: String) extends TerminalRunStatus {
    // Don't want to include errorMessage or code in the snappy status toString:
    override def toString = "Failed"
  }
}

case class Run(runId: String, pipeline: Pipeline, tag: String) {

  lazy val workflowId = pipeline.workflow.id
  lazy val call = pipeline.key.scope

  def status(): RunStatus = {
    val op = pipeline.genomicsService.operations().get(runId).execute

    if (op.getDone) {
      // If there's an error, generate a Failed status. Otherwise, we were successful!
      Option(op.getError) map { x => Failed(x.getCode, x.getMessage) } getOrElse Success
    } else if (op.hasStarted) {
      Running
    } else {
      Initializing
    }
  }

  def checkStatus(backendCall: JesBackendCall, previousStatus: Option[RunStatus]): RunStatus = {
    val currentStatus = status()

    if (!(previousStatus contains currentStatus)) {
      // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise
      // just use the state names.
      val prevStateName = previousStatus map { _.toString } getOrElse "-"
      Log.info(s"$tag: Status change from $prevStateName to $currentStatus")

      // Update the database state:
      val newBackendInfo = JesCallBackendInfo(Option(JesId(runId)), Option(JesStatus(currentStatus.toString)))
      globalDataAccess.updateExecutionBackendInfo(workflowId, CallKey(call, pipeline.key.index), newBackendInfo)

      // If this has transitioned to a running or complete state from a state this is not running or complete,
      // register the abort function.
      if (currentStatus.isRunningOrComplete && (previousStatus.isEmpty || !previousStatus.get.isRunningOrComplete)) {
        backendCall.callAbortRegistrationFunction.register(AbortFunction(() => abort()))
      }
    }
    currentStatus
  }

  def abort(): Unit = {
    val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
    pipeline.genomicsService.operations().cancel(runId, cancellationRequest).execute
  }
}