package cromwell.backend.impl.jes

import java.time.OffsetDateTime
import java.util.{ArrayList => JArrayList}

import com.google.api.client.util.{ArrayMap => GArrayMap}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.impl.jes.RunStatus.{Failed, Initializing, Running, Success}
import cromwell.core.ExecutionEvent
import cromwell.core.logging.JobLogger
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.language.postfixOps

object Run {
  private val GenomicsScopes = List(
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava

  private val AcceptableEvents = Set("start", "pulling-image", "localizing-files", "running-docker", "delocalizing-files", "ok", "fail", "start-shutdown", "preempted")

  val NoAddressFieldName = "noAddress"

  val slf4jLogger = LoggerFactory.getLogger(Run.getClass)

  def apply(runIdForResumption: Option[String],
            jobDescriptor: BackendJobDescriptor,
            runtimeAttributes: JesRuntimeAttributes,
            callRootPath: String,
            commandLine: String,
            logFileName: String,
            jesParameters: Seq[JesParameter],
            projectId: String,
            computeServiceAccount: String,
            preemptible: Boolean,
            genomicsInterface: Genomics): Run = {
    val logger = new JobLogger("JesRun", jobDescriptor.workflowDescriptor.id, jobDescriptor.key.tag, None, Set(slf4jLogger))

    lazy val workflow = jobDescriptor.workflowDescriptor
    val pipelineInfoBuilder = if (preemptible) PreemptibleJesPipelineInfoBuilder else NonPreemptibleJesPipelineInfoBuilder
    val pipelineInfo = pipelineInfoBuilder.build(commandLine, runtimeAttributes)

    val pipeline = new Pipeline()
             .setProjectId(projectId)
             .setDocker(pipelineInfo.docker)
             .setResources(pipelineInfo.resources)
             .setName(workflow.workflow.unqualifiedName)
             .setInputParameters(jesParameters.collect({ case i: JesInput => i.toGooglePipelineParameter }).toVector.asJava)
             .setOutputParameters(jesParameters.collect({ case i: JesFileOutput => i.toGooglePipelineParameter }).toVector.asJava)

    // disks cannot have mount points at runtime, so set them null
    val runtimePipelineResources = {
      val resources = pipelineInfoBuilder.build(commandLine, runtimeAttributes).resources
      val disksWithoutMountPoint = resources.getDisks.asScala map { _.setMountPoint(null) }
      resources.setDisks(disksWithoutMountPoint.asJava)
    }

    def runPipeline: String = {
      val svcAccount = new ServiceAccount().setEmail(computeServiceAccount).setScopes(GenomicsScopes)
      val rpargs = new RunPipelineArgs().setProjectId(projectId).setServiceAccount(svcAccount).setResources(runtimePipelineResources)

      rpargs.setInputs(jesParameters.collect({ case i: JesInput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
      logger.debug(s"Inputs:\n${stringifyMap(rpargs.getInputs.asScala.toMap)}")

      rpargs.setOutputs(jesParameters.collect({ case i: JesFileOutput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
      logger.debug(s"Outputs:\n${stringifyMap(rpargs.getOutputs.asScala.toMap)}")

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
    if (runIdForResumption.isDefined) { logger.info(s"JES Run is resuming with Run ID: ${runIdForResumption.get.toString}") }
    new Run(runId, genomicsInterface)
  }

  private def stringifyMap(m: Map[String, String]): String = m map { case(k, v) => s"  $k -> $v"} mkString "\n"

  implicit class RunOperationExtension(val operation: Operation) extends AnyVal {
    def hasStarted = operation.getMetadata.asScala.get("startTime") isDefined
  }

  def interpretOperationStatus(op: Operation): RunStatus = {
    if (op.getDone) {
      lazy val eventList = getEventList(op)
      lazy val ceInfo = op.getMetadata.get ("runtimeMetadata").asInstanceOf[GArrayMap[String,Object]].get("computeEngine").asInstanceOf[GArrayMap[String, String]]
      lazy val machineType = Option(ceInfo.get("machineType"))
      lazy val instanceName = Option(ceInfo.get("instanceName"))
      lazy val zone = Option(ceInfo.get("zone"))

      // If there's an error, generate a Failed status. Otherwise, we were successful!
      Option(op.getError) match {
        case None => Success(eventList, machineType, zone, instanceName)
        case Some(error) => Failed(error.getCode, Option(error.getMessage).toList, eventList, machineType, zone, instanceName)
      }
    } else if (op.hasStarted) {
      Running
    } else {
      Initializing
    }
  }

  def getEventList(op: Operation): Seq[ExecutionEvent] = {
    val metadata = op.getMetadata.asScala.toMap

    val starterEvents: Seq[ExecutionEvent] = Seq(
      eventIfExists("createTime", metadata, "waiting for quota"),
      eventIfExists("startTime", metadata, "initializing VM")).flatten

    val eventsList = for {
      events <- metadata.get("events").toSeq
      entry <- events.asInstanceOf[JArrayList[GArrayMap[String, String]]].asScala
    } yield ExecutionEvent(entry.get("description"), OffsetDateTime.parse(entry.get("startTime")))

    val filteredEventsList: Seq[ExecutionEvent] = eventsList filter { i => AcceptableEvents.contains(i.name) }

    // A little bit ugly... the endTime of the jes operation can actually be before the final "event" time, due to differences
    // in the reported precision. As a result, we have to make sure it all lines up nicely:
    val finalEvent = getCromwellPollIntervalEvent(metadata, filteredEventsList)

    starterEvents ++ filteredEventsList :+ finalEvent
  }

  private def getCromwellPollIntervalEvent(metadata: Map[String, AnyRef], eventsList: Seq[ExecutionEvent]) = {
    {
      val jesReportedEndTime = eventIfExists("endTime", metadata, "cromwell poll interval")
      val finalEventsListTime = if (eventsList.nonEmpty) Some(eventsList.last.offsetDateTime) else None

      (jesReportedEndTime, finalEventsListTime) match {
        case (Some(jesEndTime), Some(finalEventTime)) =>
          if (jesEndTime.offsetDateTime isAfter finalEventTime) jesEndTime else jesEndTime.copy(offsetDateTime = finalEventTime)
        case (Some(jesEndTime), None) => jesEndTime
        case (None, Some(finalEventTime)) => ExecutionEvent("cromwell poll interval", finalEventTime)
        case (None, None) =>
          throw new IllegalArgumentException("Both jesReportedEndTime and finalEventsListTime were None.")
      }
    }
  }

  private def eventIfExists(name: String, metadata: Map[String, AnyRef], eventName: String): Option[ExecutionEvent] = {
    metadata.get(name) map { time => ExecutionEvent(eventName, OffsetDateTime.parse(time.toString)) }
  }
}

case class Run(runId: String, genomicsInterface: Genomics) {

  def getOperationCommand = genomicsInterface.operations().get(runId)

  def abort(): Unit = {
    val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
    genomicsInterface.operations().cancel(runId, cancellationRequest).execute
    ()
  }
}
