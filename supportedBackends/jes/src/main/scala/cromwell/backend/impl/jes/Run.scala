package cromwell.backend.impl.jes

import java.time.OffsetDateTime
import java.util.{ArrayList => JArrayList}

import com.google.api.client.util.{ArrayMap => GArrayMap}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model._
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.impl.jes.RunStatus.{Failed, Initializing, Running, Success}
import cromwell.core.logging.JobLogger
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.language.postfixOps

object Run {
  private val GenomicsScopes = List(
    "https://www.googleapis.com/auth/genomics",
    "https://www.googleapis.com/auth/compute"
  ).asJava

  private val JesServiceAccount = new ServiceAccount().setEmail("default").setScopes(GenomicsScopes)
  private lazy val MaximumPollingInterval = Duration(ConfigFactory.load.getConfig("backend").getConfig("jes").getInt("maximumPollingInterval"), "seconds")
  private val InitialPollingInterval = 5 seconds
  private val PollingBackoffFactor = 1.1

  val slf4jLogger = LoggerFactory.getLogger(Run.getClass)

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
    val logger = new JobLogger("JesRun", jobDescriptor.descriptor.id, jobDescriptor.key.tag, None, Set(slf4jLogger))

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
      val eventList = getEventList(op)
      Option(op.getError) match {
        case None => Success(eventList)
        case Some(error) => Failed(error.getCode, Option(error.getMessage), eventList)
      }
    } else if (op.hasStarted) {
      Running
    } else {
      Initializing
    }
  }

  private def getEventList(op: Operation): Seq[EventStartTime] = {
    val metadata = op.getMetadata.asScala.toMap

    val starterEvents: Seq[EventStartTime] = Seq(
      eventIfExists("createTime", metadata, "waiting for quota"),
      eventIfExists("startTime", metadata, "initializing VM")).flatten

    val eventsList = for {
      events <- metadata.get("events").toSeq
      entry <- events.asInstanceOf[JArrayList[GArrayMap[String, String]]].asScala
    } yield EventStartTime(entry.get("description"), OffsetDateTime.parse(entry.get("startTime")))

    // A little bit ugly... the endTime of the jes operation can actually be before the final "event" time, due to differences
    // in the reported precision. As a result, we have to make sure it all lines up nicely:
    val finalEvent = getCromwellPollIntervalEvent(metadata, eventsList)

    starterEvents ++ eventsList :+ finalEvent
  }

  private def getCromwellPollIntervalEvent(metadata: Map[String, AnyRef], eventsList: Seq[EventStartTime]) = {
    {
      val jesReportedEndTime = eventIfExists("endTime", metadata, "cromwell poll interval")
      val finalEventsListTime = if (eventsList.nonEmpty) Some(eventsList.last.offsetDateTime) else None

      (jesReportedEndTime, finalEventsListTime) match {
        case (Some(jesEndTime), Some(finalEventTime)) =>
          if (jesEndTime.offsetDateTime isAfter finalEventTime) jesEndTime else jesEndTime.copy(offsetDateTime = finalEventTime)
        case (Some(jesEndTime), None) => jesEndTime
        case (None, Some(finalEventTime)) => EventStartTime("cromwell poll interval", finalEventTime)
        case (None, None) =>
          throw new IllegalArgumentException("Both jesReportedEndTime and finalEventsListTime were None.")
      }
    }
  }

  private def eventIfExists(name: String, metadata: Map[String, AnyRef], eventName: String): Option[EventStartTime] = {
    metadata.get(name) map {
      case time => EventStartTime(eventName, OffsetDateTime.parse(time.toString))
    }
  }

  def abort(): Unit = {
    val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
    genomicsInterface.operations().cancel(runId, cancellationRequest).execute
  }
}
