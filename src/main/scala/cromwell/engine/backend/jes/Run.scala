package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.{CancelOperationRequest, Logging, RunPipelineRequest, ServiceAccount, _}
import cromwell.engine.backend.jes.JesBackend.JesParameter
import cromwell.engine.backend.jes.Run.{Failed, Running, Success, _}
import cromwell.util.google.GoogleScopes
import org.slf4j.LoggerFactory

import scala.annotation.tailrec
import scala.collection.JavaConverters._
import scala.language.postfixOps

object Run  {
  val JesServiceAccount = new ServiceAccount().setEmail("default").setScopes(GoogleScopes.Scopes.asJava)
  lazy val Log = LoggerFactory.getLogger("main")

  def apply(pipeline: Pipeline): Run = {
    val rpr = new RunPipelineRequest().setPipelineId(pipeline.id).setProjectId(pipeline.projectId).setServiceAccount(JesServiceAccount)
    val tag = s"JES Run [UUID(${pipeline.workflow.shortId}):${pipeline.call.name}]"

    rpr.setInputs(pipeline.jesParameters.filter(_.isInput).toRunMap)
    Log.info(s"$tag Inputs:\n${stringifyMap(rpr.getInputs.asScala.toMap)}")

    rpr.setOutputs(pipeline.jesParameters.filter(_.isOutput).toRunMap)
    Log.info(s"$tag Outputs:\n${stringifyMap(rpr.getOutputs.asScala.toMap)}")

    val logging = new Logging()
    logging.setGcsPath(pipeline.gcsPath)
    rpr.setLogging(logging)

    // Currently, disk resources need to be specified both at pipeline creation and pipeline run time
    val resources = new Resources()
    resources.setDisks( scala.collection.JavaConversions.seqAsJavaList(pipeline.call.task.runtimeAttributes.defaultDisks))
    rpr.setResources(resources)

    val id = pipeline.genomicsService.pipelines().run(rpr).execute().getName
    Log.info(s"$tag JES ID is $id")
    new Run(id, pipeline, tag)
  }

  private def stringifyMap(m: Map[String, String]): String = m map { case(k, v) => s"  $k -> $v"} mkString("\n")

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
}

case class Run(name: String, pipeline: Pipeline, tag: String) {

  def status(): RunStatus = {
    val op = pipeline.genomicsService.operations().get(name).execute

    if (op.getDone) {
      // If there's an error, generate a Failed status. Otherwise, we were successful!
      Option(op.getError) map { x => Failed(x.getCode, x.getMessage) } getOrElse Success
    } else if (op.hasStarted) {
      Running
    } else {
      Initializing
    }
  }

  @tailrec
  final def waitUntilComplete(previousStatus: Option[RunStatus]): TerminalRunStatus = {
    val currentStatus = status()

    if (!(previousStatus contains currentStatus)) {
      // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise
      // just use the state names.
      val prevStateName = previousStatus map { _.toString } getOrElse "-"
      Log.info(s"$tag: Status change from $prevStateName to $currentStatus")
    }

    currentStatus match {
      case x: TerminalRunStatus => x
      case _ =>
        Thread.sleep(5000)
        waitUntilComplete(Option(currentStatus))
    }
  }

  @tailrec
  final def waitUntilRunningOrComplete(): Unit = {
    val currentStatus = status()
    currentStatus match {
      case Initializing => // Done
      case x: TerminalRunStatus => // Done
      case _ =>
        Thread.sleep(5000)
        waitUntilRunningOrComplete()
    }
  }

  def abort(): Unit = {
    val cancellationRequest: CancelOperationRequest = new CancelOperationRequest()
    pipeline.genomicsService.operations().cancel(name, cancellationRequest).execute
  }
}
