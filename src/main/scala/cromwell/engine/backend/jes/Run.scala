package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.{ServiceAccount, RunPipelineRequest}
import cromwell.engine.backend.jes.JesBackend.JesParameter
import cromwell.engine.backend.jes.Run.{Running, Success, Failed}
import cromwell.util.google.GoogleScopes
import scala.annotation.tailrec
import scala.collection.JavaConverters._

import Run._

object Run {
  val JesServiceAccount = new ServiceAccount().setEmail("default").setScopes(GoogleScopes.Scopes.asJava)

  def apply(pipeline: Pipeline): Run = {
    val rpr = new RunPipelineRequest().setPipelineId(pipeline.id).setProjectId(pipeline.projectId).setServiceAccount(JesServiceAccount)

    rpr.setInputs(pipeline.jesParameters.filter(_.isInput).toRunMap)
    println(s"Run inputs are ${rpr.getInputs}")
    rpr.setOutputs(pipeline.jesParameters.filter(_.isOutput).toRunMap)
    println(s"Run outputs are ${rpr.getOutputs}")

    val id = pipeline.genomicsService.pipelines().run(rpr).execute().getName
    println(s"Run Id is $id")
    new Run(id, pipeline)
  }

  implicit class RunJesParameters(val params: Seq[JesParameter]) extends AnyVal {
    def toRunMap = params.map(p => p.name -> p.gcs).toMap.asJava
  }

  sealed trait RunStatus // FIXME: These dates shouldn't be Strings
  trait TerminalRunStatus extends RunStatus
  final case class Initializing(created: String) extends RunStatus
  final case class Running(created: String, started: String) extends RunStatus
  final case class Success(created: String, started: String, finished: String) extends TerminalRunStatus
  // FIXME: Not capturing the errorDetails map, might want to do that (but .asScala issues exist there))
  // FIXME: Also it's been a while since I looked but there never seemed to be anything in there anyways
  final case class Failed(created: String, started: String, finished: String, errorCode: Int, errorMessage: String) extends TerminalRunStatus
}

case class Run(name: String, pipeline: Pipeline) {
  def status(): RunStatus = {
    val op = pipeline.genomicsService.operations().get(name).execute
    val metadata = OperationMetadata(op)

    if (op.getDone) {
      val error = Option(op.getError)
      error match { // "What is this, amateur hour? ...."
        case Some(x) => Failed(
          metadata.created,
          /* metadata.started.get*/ "started",
          /* metadata.finished.get */ "finished",
          x.getCode,
          x.getMessage)
        case None => Success(metadata.created, metadata.started.get, metadata.finished.get)
      }
    } else metadata.started map {s => Running(metadata.created, s)} getOrElse Initializing(metadata.created)
  }

  @tailrec
  final def waitUntilComplete(): TerminalRunStatus = {
    val currentStatus = status()
    println(s"Current status is $currentStatus")
    currentStatus match {
      case x: TerminalRunStatus => x
      case _ =>
        Thread.sleep(5000)
        waitUntilComplete()
    }
  }
}