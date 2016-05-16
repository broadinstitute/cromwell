package centaur.test.standard

import java.nio.file.{Path, Paths}

import cats.data.Validated._
import centaur.test.ErrorOr
import com.typesafe.config.Config
import lenthall.config.ScalaConfig._

// TODO: metadata will disappear in a few PRs
case class StandardTestCaseFiles(wdl: Path, inputs: Option[Path], options: Option[Path], metadata: Option[Path])

object StandardTestCaseFiles {
  def fromConfig(conf: Config, basePath: Path): ErrorOr[StandardTestCaseFiles] = {
    conf.getPathOption("wdl") match {
      case Some(wdl) => Valid(StandardTestCaseFiles(basePath.resolve(wdl), conf, basePath))
      case None => invalidNel("No wdl path provided")
    }
  }

  def apply(wdl: Path, conf: Config, basePath: Path): StandardTestCaseFiles = {
    def getOptionalPath(name: String) = conf.getPathOption(name) map basePath.resolve
    StandardTestCaseFiles(wdl, getOptionalPath("inputs"), getOptionalPath("options"), getOptionalPath("metadata"))
  }

  implicit class EnhancedConfig(val conf: Config) extends AnyVal {
    def getPathOption(key: String): Option[Path] = conf.getStringOption(key) map { Paths.get(_) }
  }
}
