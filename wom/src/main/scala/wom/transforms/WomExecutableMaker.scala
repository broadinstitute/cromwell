package wom.transforms

import common.Checked
import wom.executable.Executable
import simulacrum._
import wom.callable.{TaskDefinition, WorkflowDefinition}
import wom.transforms.WomExecutableMaker.ExecutableMakerInputs

import scala.concurrent.Future
import scala.language.implicitConversions

@typeclass
trait WomExecutableMaker[A] {
  def toWomExecutable(inputs: ExecutableMakerInputs[A]): Checked[Executable]
  def toWomExecutable(a: A, inputs: Option[String] = None): Checked[Executable]
}

object WomExecutableMaker {
  final case class ImportResolution(workfows: WorkflowDefinition, tasks: TaskDefinition)
  final case class ExecutableMakerInputs[A](from: A, importResolvers: List[String => Future[Checked[ImportResolution]]], inputs: Option[String])
}
