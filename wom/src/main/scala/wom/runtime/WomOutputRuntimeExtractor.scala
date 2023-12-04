package wom.runtime

import cats.syntax.apply._
import com.typesafe.config.Config
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

object WomOutputRuntimeExtractor {
  def fromConfig(config: Config) = {
    val dockerImage = validate(config.getAs[String]("docker-image"))
    val command = validate(config.as[String]("command"))

    (dockerImage, command).mapN(WomOutputRuntimeExtractor.apply)
  }
}

case class WomOutputRuntimeExtractor(dockerImage: Option[String], command: String)
