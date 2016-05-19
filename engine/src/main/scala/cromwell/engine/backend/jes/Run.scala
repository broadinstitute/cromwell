package cromwell.engine.backend.jes

import com.google.api.client.util.ArrayMap
import com.google.api.services.genomics.{Genomics, model}
import com.google.api.services.genomics.model.{CancelOperationRequest, LoggingOptions, Pipeline, RunPipelineArgs, RunPipelineRequest, ServiceAccount, _}
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import cromwell.engine.backend.BackendCallJobDescriptor
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.jes.Run.{Failed, Running, Success, _}
import cromwell.engine.db.DataAccess._
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{AbortFunction, ExecutionEventEntry}
import cromwell.logging.WorkflowLogger
import cromwell.util.google.GoogleScopes
import org.joda.time.DateTime
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

object Run  {
  val JesServiceAccount = new ServiceAccount().setEmail("default").setScopes(GoogleScopes.Scopes.asJava)
  lazy val MaximumPollingInterval = Duration(ConfigFactory.load.getConfig("backend").getConfig("jes").getInt("maximumPollingInterval"), "seconds")
  val InitialPollingInterval = 5 seconds
  val PollingBackoffFactor = 1.1

  def apply(runIdForResumption: Option[String],
            jesJobDescriptor: JesJobDescriptor,
            jesParameters: Seq[JesParameter],
            projectId: String,
            genomicsInterface: Genomics): Run = {
    lazy val jobDescriptor = jesJobDescriptor.jobDescriptor
    lazy val workflow = jobDescriptor.workflowDescriptor
    lazy val command = jesJobDescriptor.jesCommandLine
    lazy val runtimeAttributes = jobDescriptor.callRuntimeAttributes
    lazy val key = jobDescriptor.key
    lazy val gcsPath = jobDescriptor.callRootPath.toString

    val logger = WorkflowLogger(
      "JES Run",
      workflow,
      otherLoggers = Seq(LoggerFactory.getLogger(getClass.getName)),
      callTag = Option(key.tag)
    )

    logger.debug(s"Command line is: $command")

    val runtimeInfo = if (jesJobDescriptor.preemptible) PreemptibleJesRuntimeInfo(command, runtimeAttributes) else NonPreemptibleJesRuntimeInfo(command, runtimeAttributes)
    val pipeline = new model.Pipeline()
                    .setProjectId(projectId)
                    .setDocker(runtimeInfo.docker)
                    .setResources(runtimeInfo.resources)
                    .setName(workflow.name)
                    .setInputParameters(jesParameters.collect({ case i: JesInput => i.toGooglePipelineParameter }).toVector.asJava)
                    .setOutputParameters(jesParameters.collect({ case i: JesFileOutput => i.toGooglePipelineParameter }).toVector.asJava)

    def runPipeline(): String = {
      val rpargs = new RunPipelineArgs().setProjectId(projectId).setServiceAccount(JesServiceAccount)

      rpargs.setInputs(jesParameters.collect({ case i: JesInput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
      logger.info(s"Inputs:\n${stringifyMap(rpargs.getInputs.asScala.toMap)}")

      rpargs.setOutputs(jesParameters.collect({ case i: JesFileOutput => i.name -> i.toGoogleRunParameter }).toMap.asJava)
      logger.info(s"Outputs:\n${stringifyMap(rpargs.getOutputs.asScala.toMap)}")

      val rpr = new RunPipelineRequest().setEphemeralPipeline(pipeline).setPipelineArgs(rpargs)

      val logging = new LoggingOptions()
      logging.setGcsPath(s"$gcsPath/${JesBackend.jesLogFilename(key)}")
      rpargs.setLogging(logging)

      val runId = genomicsInterface.pipelines().run(rpr).execute().getName
      logger.info(s"JES Run ID is $runId")
      runId
    }

    // If runIdForResumption is defined use that, otherwise we'll create a new Run with an ephemeral pipeline.
    // The Run code takes care of polling.
    val runId = runIdForResumption getOrElse runPipeline
    new Run(runId, workflow.id, key, genomicsInterface, logger)
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
  case class Success(events: Seq[ExecutionEventEntry]) extends TerminalRunStatus {
    override def toString = "Success"
  }
  final case class Failed(errorCode: Int, errorMessage: Option[String], events: Seq[ExecutionEventEntry]) extends TerminalRunStatus {
    // Don't want to include errorMessage or code in the snappy status toString:
    override def toString = "Failed"
  }

  // An event with a startTime timestamp
  private case class EventStartTime(name: String, timestamp: DateTime)

  def getEventList(op: Operation): Seq[ExecutionEventEntry] = {
    val starterEvents = eventIfExists("createTime", op, "waiting for quota") ++ eventIfExists("startTime", op, "initializing VM")

    val eventsList: Seq[EventStartTime] = if (op.getMetadata.containsKey("events")) {
      op.getMetadata.get("events").asInstanceOf[java.util.ArrayList[AnyRef]].asScala map { x =>
        val entry = x.asInstanceOf[ArrayMap[String, String]]
        EventStartTime(entry.get("description"), DateTime.parse(entry.get("startTime")))
      } toSeq
    } else Seq.empty

    // The final event is only used as the book-end for the final pairing (see below) so the name is never actually used...
    // ... which is rather a pity actually - it's a jolly good name.
    val finaleEvents = eventIfExists("endTime", op, "cromwell poll interval") ++ Seq(
      EventStartTime("The Queen flying around with a jet-pack, with Winston Churchill cheering and waving a huge Union Jack in the background", DateTime.now))

    // Join the Seqs together, pair up consecutive elements then make events with start and end times.
    ((starterEvents ++ eventsList ++ finaleEvents).sliding(2) toSeq) map { case Seq(a, b) => ExecutionEventEntry(a.name, a.timestamp, b.timestamp) }
  }

  private def eventIfExists(name: String, op: Operation, eventName: String): Seq[EventStartTime] = {
    val metadata = op.getMetadata
    if(metadata.containsKey(name))
      Seq(EventStartTime(eventName, DateTime.parse(metadata.get(name).asInstanceOf[String])))
    else
      Seq.empty
  }
}

case class Run(runId: String, workflowId: WorkflowId, key: BackendCallKey, genomicsInterface: Genomics, logger: WorkflowLogger) {
  lazy val call = key.scope

  def status(): RunStatus = {
    val op = genomicsInterface.operations().get(runId).execute
    if (op.getDone) {
      // If there's an error, generate a Failed status. Otherwise, we were successful!
      val eventList = getEventList(op)
      Option(op.getError) map { x => Failed(x.getCode, Option(x.getMessage), eventList) } getOrElse Success(eventList)
    } else if (op.hasStarted) {
      Running
    } else {
      Initializing
    }
  }

  def checkStatus(jobDescriptor: BackendCallJobDescriptor, previousStatus: Option[RunStatus]): RunStatus = {
    val currentStatus = status()

    if (!(previousStatus contains currentStatus)) {
      // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise
      // just use the state names.
      val prevStateName = previousStatus map { _.toString } getOrElse "-"
      logger.info(s"Status change from $prevStateName to $currentStatus")

      /*
      TODO: Not sure we're supposed to be directly talking to the database.
      This doesn't even wait for the future to complete. Pretty sure this should be a message to the workflow actor,
      that then contacts the database to change the state. For now, updating this end run to the database to pass in the
      default, global execution context.
       */

      // Update the database state:
      // TODO the database API should probably be returning DBIOs so callers can compose and wrap with a transaction.
      globalDataAccess.updateExecutionInfo(workflowId, BackendCallKey(call, key.index, key.attempt), JesBackend.InfoKeys.JesRunId, Option(runId))(ExecutionContext.global)
      globalDataAccess.updateExecutionInfo(workflowId, BackendCallKey(call, key.index, key.attempt), JesBackend.InfoKeys.JesStatus, Option(currentStatus.toString))(ExecutionContext.global)

      // If this has transitioned to a running or complete state from a state that is not running or complete,
      // register the abort function.
      if (currentStatus.isRunningOrComplete && (previousStatus.isEmpty || !previousStatus.get.isRunningOrComplete)) {
        jobDescriptor.abortRegistrationFunction.foreach(_.register(AbortFunction(() => abort())))
      }
    }

    currentStatus
  }

  def abort(): Unit = {
    val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
    genomicsInterface.operations().cancel(runId, cancellationRequest).execute
  }
}
