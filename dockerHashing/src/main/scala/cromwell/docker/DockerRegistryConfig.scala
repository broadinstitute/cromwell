package cromwell.docker

import java.util.concurrent.Executors

import cats.syntax.apply._
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.core.io.Throttle

import scala.concurrent.ExecutionContext

case class DockerRegistryConfig(throttle: Option[Throttle], nbThreads: Int) {
  lazy val executionContext = ExecutionContext.fromExecutor(Executors.newFixedThreadPool(nbThreads))
}

object DockerRegistryConfig {
  lazy val default = DockerRegistryConfig(None, 5)

  def fromConfig(config: Config): ErrorOr[DockerRegistryConfig] = {
    val throttle = validate { config.getAs[Throttle]("throttle") }
    val threads = validate { config.as[Int]("num-threads") }
   
    (throttle, threads) mapN DockerRegistryConfig.apply
  }
}
