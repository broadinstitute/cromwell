package centaur.test

import cats.data.Validated._
import com.typesafe.config.Config
import configs.Result.{Failure, Success}
import configs.syntax._

object TagsStrings {

  def fromConfig(conf: Config): ErrorOr[List[String]] = {
    conf.get[List[String]]("tags") match {
      case Success(tagStrings) => Valid(tagStrings.map(_.toLowerCase).distinct)
      case Failure(_) => Valid(List.empty[String])
    }
  }
}
