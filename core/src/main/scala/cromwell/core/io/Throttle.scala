package cromwell.core.io

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ValueReader

import scala.concurrent.duration.FiniteDuration

case class Throttle(elements: Int, per: FiniteDuration, maximumBurst: Int) {
  lazy val delay: FiniteDuration = per.div(elements.toLong)
}

object Throttle {
  implicit val throttleOptionValueReader: ValueReader[Option[Throttle]] = (config: Config, path: String) => {
    config.getAs[Config](path) map { throttleConfig =>
      val elements = throttleConfig.as[Int]("number-of-requests")
      val per = throttleConfig.as[FiniteDuration]("per")
      Throttle(elements, per, elements)
    }
  }
}
