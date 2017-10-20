package centaur.test

import configs.Result.{Failure, Success}
import configs.syntax._
import cats.data.Validated._
import cats.Apply
import com.typesafe.config.Config
import lenthall.validation.ErrorOr.ErrorOr

case class TestOptions(tags: List[String], ignore: Boolean)

object TestOptions {

  def fromConfig(conf: Config): ErrorOr[TestOptions] = {
    val tags = tagsFromConfig(conf)
    val ignore = ignoreFromConfig(conf)

    Apply[ErrorOr].map2(tags, ignore)((t, i) => TestOptions(t, i))
  }

  def tagsFromConfig(conf: Config): ErrorOr[List[String]] = {
    conf.get[List[String]]("tags") match {
      case Success(tagStrings) => Valid(tagStrings.map(_.toLowerCase).distinct)
      case Failure(_) => Valid(List.empty[String])
    }
  }

  def ignoreFromConfig(conf: Config): ErrorOr[Boolean] = {

    if (conf.hasPath("ignore")) {
      conf.get[Boolean]("ignore") match {
        case Success(ignore) => Valid(ignore)
        case Failure(f) => invalidNel(s"Invalid 'ignore' value: $f")
      }
    } else {
      Valid(false)
    }
  }
}

