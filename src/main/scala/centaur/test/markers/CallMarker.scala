package centaur.test.markers

import cats.syntax.validated._
import com.typesafe.config.Config
import configs.Result.{Failure, Success}
import configs.syntax._
import lenthall.validation.ErrorOr.ErrorOr

object CallMarker {
  def fromConfig(config: Config): ErrorOr[Option[CallMarker]] = {
    config.get[Option[String]]("callMark") match {
      case Success(marker) => (marker map CallMarker.apply).validNel
      case Failure(f) => s"Invalid restart marker $f".invalidNel
    }
  }
}

/**
  * Allows to mark a call in a test description for which a specific action should be performed
  * Test formulas need to explicitly handle those markers in order for them to be honored
  */
case class CallMarker(callKey: String)
