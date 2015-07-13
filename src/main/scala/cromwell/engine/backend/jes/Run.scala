package cromwell.engine.backend.jes

import com.google.api.services.genomics.model.RunPipelineRequest
import cromwell.engine.backend.jes.Run.{Running, Success, Failed}
import scala.annotation.tailrec
import scala.collection.JavaConverters._

import Run._

object Run {
  def apply(pipeline: Pipeline): Run = {
    val rpr = new RunPipelineRequest().setPipelineId(pipeline.id).setProjectId(pipeline.projectId).setServiceAccount(JesBackend.JesServiceAccount)

    /*
       FIXME: This is seriously crappy and almost certainly isn't the right way to even go about things going forward even if it was pretty

       We have the root filename (e.g. 'foo.txt') from inputsToLocalize, so we can slap on this call's GCS path.
       However that file came from a different call - but we don't directly know from where. However we have access
       to this call and can see it's input mappings, and the two share the same key. So this is grabbing the original
       call name from those input mappings and replacing the current call name for that.
    */
    rpr.setInputs(pipeline.inputsToLocalize.map {case (k, v) =>
      val blah = s"${pipeline.gcsPath}/${v.toString}"
      val origCall = pipeline.call.inputMappings(k).toString.split("\\.").head
      val filePath = blah.replace(pipeline.call.name, origCall)
      k -> filePath
    }.asJava)
    println(s"Run inputs are ${rpr.getInputs}")

    // FIXME: Outputs - these are currently hardcoded, obviously. We'll also need the JesEngineFunctions before outputs are particularly useful
    val outputs = Map("stdout" -> s"${pipeline.gcsPath}/stdout.txt", "stderr" -> s"${pipeline.gcsPath}/stderr.txt")
    println(s"Run outputs are $outputs")
    rpr.setOutputs(outputs.asJava)

    val id = pipeline.genomicsService.pipelines().run(rpr).execute().getName
    println(s"Run Id is $id")
    new Run(id, pipeline)
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
        case Some(x) => Failed(metadata.created, metadata.started.get, metadata.finished.get, x.getCode, x.getMessage)
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