package cromwell.core.path

import java.nio.file.Path

import scala.util.Try

trait PathBuilder {
  def name: String
  def build(pathAsString: String): Try[Path]
}
