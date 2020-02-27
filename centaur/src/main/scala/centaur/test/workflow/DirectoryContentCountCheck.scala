package centaur.test.workflow

import cats.data.Validated._
import cats.implicits._
import centaur.test.{AWSFilesChecker, FilesChecker, LocalFilesChecker, PipelinesFilesChecker}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import configs.Result
import configs.syntax._

final case class DirectoryContentCountCheck(expectedDirectoryContentsCounts: Map[String, Int], checkFiles: FilesChecker)

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

      val fileSystemChecker: ErrorOr[FilesChecker] = conf.get[String]("fileSystemCheck") match {
        case Result.Success("gcs") => valid(PipelinesFilesChecker)
        case Result.Success("local") => valid(LocalFilesChecker)
        case Result.Success("aws") => valid(AWSFilesChecker)
        case Result.Success(_) => invalidNel(s"Test '$name': Invalid 'fileSystemCheck' value (must be either 'local', 'gcs' or 'aws'")
        case Result.Failure(_) => invalidNel(s"Test '$name': Must specify a 'fileSystemCheck' value (must be either 'local', 'gcs' or 'aws'")
      }

      (directoryContentCountsValidation, fileSystemChecker) mapN { (d, f) => Option(DirectoryContentCountCheck(d, f)) }
    }
  }
}
