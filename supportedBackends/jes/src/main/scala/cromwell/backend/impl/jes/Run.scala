package cromwell.backend.impl.jes

import com.google.api.services.genomics.model._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobDescriptor
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.language.postfixOps

object Run  {
  private val GenomicsScopes = List(
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava

  private val JesServiceAccount = new ServiceAccount().setEmail("default").setScopes(GenomicsScopes)
  private lazy val MaximumPollingInterval = Duration(ConfigFactory.load.getConfig("backend").getConfig("jes").getInt("maximumPollingInterval"), "seconds")
  private val InitialPollingInterval = 5 seconds
  private val PollingBackoffFactor = 1.1

  def apply(pipeline: Pipeline, logFileName: String): Run = {
    val logger = LoggerFactory.getLogger(Run.getClass)

    if (pipeline.pipelineId.isDefined == pipeline.runIdForResumption.isDefined) {
      val message =
        s"""
          |Exactly one of JES pipeline ID or run ID for resumption must be specified to create a Run.
          |pipelineId = ${pipeline.pipelineId}, runIdForResumption = ${pipeline.runIdForResumption}.
        """.stripMargin
      throw new RuntimeException(message)
    }

    def runPipeline: String = {
      val rpargs = new RunPipelineArgs().setProjectId(pipeline.projectId).setServiceAccount(JesServiceAccount)

      rpargs.setInputs(pipeline.jesParameters.collect({ case i: JesInput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
      logger.info(s"Inputs:\n${stringifyMap(rpargs.getInputs.asScala.toMap)}")

      rpargs.setOutputs(pipeline.jesParameters.collect({ case i: JesFileOutput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
      logger.info(s"Outputs:\n${stringifyMap(rpargs.getOutputs.asScala.toMap)}")

      val rpr = new RunPipelineRequest().setPipelineId(pipeline.pipelineId.get).setPipelineArgs(rpargs)

      val logging = new LoggingOptions()
      logging.setGcsPath(s"${pipeline.gcsPath}/$logFileName")
      rpargs.setLogging(logging)

      val runId = pipeline.genomicsService.pipelines().run(rpr).execute().getName
      logger.info(s"JES Run ID is $runId")
      runId
    }

    // Only run the pipeline if the pipeline ID is defined.  The pipeline ID not being defined corresponds to a
    // resumption of a previous run, and runIdForResumption will be defined.  The Run code takes care of polling
    // in both the newly created and resumed scenarios.
    val runId = if (pipeline.pipelineId.isDefined) runPipeline else pipeline.runIdForResumption.get
    new Run(runId, pipeline /*, logger */)
  }

  private def stringifyMap(m: Map[String, String]): String = m map { case(k, v) => s"  $k -> $v"} mkString "\n"

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
  case object Success extends TerminalRunStatus {
    override def toString = "Success"
  }
  final case class Failed(errorCode: Int, errorMessage: Option[String]) extends TerminalRunStatus {
    // Don't want to include errorMessage or code in the snappy status toString:
    override def toString = "Failed"
  }

  // PBE deleted all the ExecutionEvent stuff
}

// PBE hacked out loggers
case class Run(runId: String, pipeline: Pipeline /*, logger: WorkflowLogger */) {
  import Run._

  lazy val workflowId = pipeline.jobDescriptor.descriptor.id
  lazy val call = pipeline.jobDescriptor.key.scope

  def status(): RunStatus = {
    val op = pipeline.genomicsService.operations().get(runId).execute
    if (op.getDone) {
      // If there's an error, generate a Failed status. Otherwise, we were successful!
      // val eventList = getEventList(op)
      Option(op.getError) map { x => Failed(x.getCode, Option(x.getMessage) /*, eventList */) } getOrElse Success
    } else if (op.hasStarted) {
      Running
    } else {
      Initializing
    }
  }

  def checkStatus(jobDescriptor: BackendJobDescriptor, previousStatus: Option[RunStatus]): RunStatus = {
    val currentStatus = status()

    if (!(previousStatus contains currentStatus)) {
      // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise
      // just use the state names.
      val prevStateName = previousStatus map { _.toString } getOrElse "-"
      println(s"Status change from $prevStateName to $currentStatus")

      // PBE deleted db writes for run id and status
      // PBE deleted abort function registration
    }
    currentStatus
  }

  def abort(): Unit = {
    val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
    pipeline.genomicsService.operations().cancel(runId, cancellationRequest).execute
  }
}
