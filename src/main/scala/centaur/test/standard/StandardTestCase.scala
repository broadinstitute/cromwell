package centaur.test.standard

import java.nio.file.Path

import cats.Apply
import cats.data.Validated._
import cats.std.list._
import centaur.test._
import com.typesafe.config.{Config, ConfigFactory}
import lenthall.config.ScalaConfig._
import StandardTestCaseFiles.EnhancedConfig

import scala.util.{Failure, Success, Try}

case class StandardTestCase(name: String, basePath: Path, testFormat: StandardTestFormat, files: StandardTestCaseFiles)

object StandardTestCase {
  def fromPath(path: Path): ErrorOr[StandardTestCase] = {
    Try(ConfigFactory.parseFile(path.toFile)) match {
      case Success(c) => StandardTestCase.fromConfig(c, path.getParent)
      case Failure(f) => invalidNel(s"Invalid test config: $path")
    }
  }

  def fromConfig(conf: Config, configPath: Path): ErrorOr[StandardTestCase] = {
    conf.getStringOption("testName") match {
      case Some(n) =>
        // If basePath is provided it'll be used as basis for finding other files, otherwise use the dir the config was in
        val basePath = conf.getPathOption("basePath") getOrElse configPath
        val format = StandardTestFormat.fromConfig(conf)
        val files = conf.getConfigOption("files") match {
          case Some(c) => StandardTestCaseFiles.fromConfig(c, basePath)
          case None => invalidNel(s"No 'files' block in $configPath")
        }

        Apply[ErrorOr].map2(format, files)((a, b) => StandardTestCase(n, basePath, a, b))
      case None => invalidNel(s"No testName for: $configPath")
    }
  }
}

