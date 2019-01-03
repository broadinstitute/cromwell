package centaur.test.workflow

import cats.data.Validated._
import cats.implicits._
import centaur.test.{CheckFiles, PipelinesCheckFiles, LocalCheckFiles}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import configs.Result
import configs.syntax._

final case class DirectoryContentCountCheck(expectedDrectoryContentsCounts: Map[String, Int], checkFiles: CheckFiles)

object DirectoryContentCountCheck {
  def forConfig(name: String, conf: Config): ErrorOr[Option[DirectoryContentCountCheck]] = {
    val keepGoing = conf.hasPath("outputExpectations")
    if (!keepGoing) {
      valid(None)
    } else {
      val directoryContentCountsValidation: ErrorOr[Map[String, Int]] = conf.get[Map[String, Int]]("outputExpectations") match {
        case Result.Success(a) => valid(a)
        case Result.Failure(_) => invalidNel(s"Test '$name': Unable to read outputExpectations as a Map[String, Int]")
      }

      val fileSystemChecker: ErrorOr[CheckFiles] = conf.get[String]("fileSystemCheck") match {
        case Result.Success("gcs") => valid(PipelinesCheckFiles)
        case Result.Success("local") => valid(LocalCheckFiles)
        case Result.Success(_) => invalidNel(s"Test '$name': Invalid 'fileSystemCheck' value (must be either 'local' or 'gcs')")
        case Result.Failure(_) => invalidNel(s"Test '$name': Must specify a 'fileSystemCheck' value (must be either 'local' or 'gcs')")
      }

      (directoryContentCountsValidation, fileSystemChecker) mapN { (d, f) => Option(DirectoryContentCountCheck(d, f)) }
    }
  }
}
