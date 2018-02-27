package wom.transforms

import common.Checked
import wom.executable.{Executable, WomBundle}
import simulacrum._
import wom.core.WorkflowJson

import scala.concurrent.Future
import scala.language.implicitConversions

@typeclass
trait WomExecutableMaker[A] {
  def toWomExecutable(a: A, inputs: Option[WorkflowJson]): Checked[Executable]
}

@typeclass
trait WomBundleMaker[A] {
  def toWomBundle(a: A, importResolvers: List[String => Future[Checked[WomBundle]]]): Checked[WomBundle]
}

