package cromwell.backend.impl.jes

import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.impl.jes.RunStatus.{Failed, Initializing, Running, Success}
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

  def apply(runIdForResumption: Option[String],
            jobDescriptor: BackendJobDescriptor,
            runtimeAttributes: JesRuntimeAttributes,
            callRootPath: String,
            commandLine: String,
            logFileName: String,
            jesParameters: Seq[JesParameter],
            projectId: String,
            preemptible: Boolean,
            genomicsInterface: Genomics): Run = {
    val logger = LoggerFactory.getLogger(Run.getClass)

    lazy val workflow = jobDescriptor.descriptor
    val runtimeInfoBuilder = if (preemptible) PreemptibleJesRuntimeInfoBuilder else NonPreemptibleJesRuntimeInfoBuilder
    val runtimeInfo = runtimeInfoBuilder.build(commandLine, runtimeAttributes)

    val pipeline = new Pipeline()
             .setProjectId(projectId)
             .setDocker(runtimeInfo.docker)
             .setResources(runtimeInfo.resources)
             .setName(workflow.workflowNamespace.workflow.unqualifiedName)
             .setInputParameters(jesParameters.collect({ case i: JesInput => i.toGooglePipelineParameter }).toVector.asJava)
             .setOutputParameters(jesParameters.collect({ case i: JesFileOutput => i.toGooglePipelineParameter }).toVector.asJava)

    def runPipeline: String = {
      val rpargs = new RunPipelineArgs().setProjectId(projectId).setServiceAccount(JesServiceAccount)

      rpargs.setInputs(jesParameters.collect({ case i: JesInput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
      logger.info(s"Inputs:\n${stringifyMap(rpargs.getInputs.asScala.toMap)}")

      rpargs.setOutputs(jesParameters.collect({ case i: JesFileOutput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
      logger.info(s"Outputs:\n${stringifyMap(rpargs.getOutputs.asScala.toMap)}")

      val rpr = new RunPipelineRequest().setEphemeralPipeline(pipeline).setPipelineArgs(rpargs)

      val logging = new LoggingOptions()
      logging.setGcsPath(s"$callRootPath/$logFileName")
      rpargs.setLogging(logging)

      val runId = genomicsInterface.pipelines().run(rpr).execute().getName
      logger.info(s"JES Run ID is $runId")
      runId
    }

    // If runIdForResumption is defined use that, otherwise we'll create a new Run with an ephemeral pipeline.
    val runId = runIdForResumption getOrElse runPipeline
    new Run(runId, jobDescriptor, genomicsInterface)
  }

  private def stringifyMap(m: Map[String, String]): String = m map { case(k, v) => s"  $k -> $v"} mkString "\n"

  implicit class RunOperationExtension(val operation: Operation) extends AnyVal {
    def hasStarted = operation.getMetadata.asScala.get("startTime") isDefined
  }

  // PBE deleted all the ExecutionEvent stuff
}

// PBE hacked out loggers
case class Run(runId: String,  jobDescriptor: BackendJobDescriptor, genomicsInterface: Genomics/*, logger: WorkflowLogger */) {
  import Run._

  lazy val workflowId = jobDescriptor.descriptor.id
  lazy val call = jobDescriptor.key.scope

  def status(): RunStatus = {
    val op = genomicsInterface.operations().get(runId).execute
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
    genomicsInterface.operations().cancel(runId, cancellationRequest).execute
  }
}
