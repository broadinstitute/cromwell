package cromwell.engine.workflow.workflowstore

import java.util.UUID

import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigFactory}
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.engine.workflow.workflowstore.WorkflowHeartbeatConfig._
import io.circe._
import io.circe.generic.semiauto._
import io.circe.syntax._
import mouse.all._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

/**
  *
  * @param cromwellId         Identifier for this Cromwell for horizontaling purposes.
  * @param heartbeatInterval  How long between `WorkflowActor` heartbeats and `WorkflowStoreHeartbeatWriteActor` database flushes.
  * @param ttl                Entries in the workflow store with heartbeats older than `ttl` ago are presumed abandoned and up for grabs.
  * @param failureShutdownDuration How long writing heartbeats are allow to fail before cromwell shuts down.
  * @param writeBatchSize     The maximum size of a write batch.
  * @param writeThreshold     The threshold of heartbeat writes above which load is considered high.
  */
case class WorkflowHeartbeatConfig(
                                    cromwellId: String,
                                    heartbeatInterval: FiniteDuration,
                                    ttl: FiniteDuration,
                                    failureShutdownDuration: FiniteDuration,
                                    writeBatchSize: Int,
                                    writeThreshold: Int) {
  override def toString: String = this.asInstanceOf[WorkflowHeartbeatConfig].asJson.spaces2
}

object WorkflowHeartbeatConfig {

  // NOTE: If these are made fully private the compiler wrongly emits a warning about them being unused, which our
  // compiler flag settings then promote to an error.

  // NOTE: This is a different encoding than circe's finiteDurationEncoder: https://github.com/circe/circe/pull/978
  private[engine] implicit lazy val encodeFiniteDuration: Encoder[FiniteDuration] = {
    Encoder.encodeString.contramap(_.toString)
  }
  private[engine] implicit lazy val encodeWorkflowHeartbeatConfig: Encoder[WorkflowHeartbeatConfig] = deriveEncoder

  def apply(config: Config): WorkflowHeartbeatConfig = {
    validate(config).toTry("Errors parsing WorkflowHeartbeatConfig").get
  }

  private def validate(config: Config): ErrorOr[WorkflowHeartbeatConfig] = {
    val randomSuffix = config
      .getOrElse("system.cromwell_id_random_suffix", true)
      .fold("-" + UUID.randomUUID().toString.take(7), "")
    val cromwellId: String = config.getOrElse("system.cromwell_id", "cromid") + randomSuffix

    val heartbeats: Config = config.getOrElse("system.workflow-heartbeats", ConfigFactory.empty())
    // Use defaults from reference.conf but don't allow values less than 10 secs except for testing purposes.
    val minHeartbeatTtl = config.getOrElse("danger.debug.only.minimum-heartbeat-ttl", 10.seconds)
    val ttl: FiniteDuration = heartbeats.as[FiniteDuration]("ttl").max(minHeartbeatTtl)
    // Use the defaults in reference.conf but don't allow values less than one third of the minimum TTL.
    val heartbeatIntervalValidation: ErrorOr[FiniteDuration] = {
      val minHeartbeatInterval = minHeartbeatTtl / 3
      val heartbeatInterval = heartbeats.as[FiniteDuration]("heartbeat-interval") max minHeartbeatInterval
      if (heartbeatInterval >= ttl) {
        val errorMessage =
          s"The system.workflow-heartbeats.heartbeat-interval ($heartbeatInterval)" +
            s" is not less than the system.workflow-heartbeats.ttl ($ttl)."
        errorMessage.invalidNel
      } else {
        heartbeatInterval.valid
      }
    }

    val failureShutdownDurationValidation: ErrorOr[FiniteDuration] = {
      // Don't allow shutdown duration to be more than the heartbeat TTL.
      val failureShutdownDuration = heartbeats.as[FiniteDuration]("write-failure-shutdown-duration") max 0.seconds
      if (failureShutdownDuration > ttl) {
        val errorMessage =
          s"The system.workflow-heartbeats.write-failure-shutdown-duration ($failureShutdownDuration)" +
            s" is greater than the system.workflow-heartbeats.ttl ($ttl)."
        errorMessage.invalidNel
      } else {
        failureShutdownDuration.valid
      }
    }

    // Sanity bound all other parameters as well.
    val writeBatchSize: Int = heartbeats.getInt("write-batch-size") max 1
    val writeThreshold: Int = heartbeats.getInt("write-threshold") max 1

    (heartbeatIntervalValidation, failureShutdownDurationValidation) mapN {
      (heartbeatInterval, failureShutdownDuration) =>
        WorkflowHeartbeatConfig(
          cromwellId = cromwellId,
          heartbeatInterval = heartbeatInterval,
          ttl = ttl,
          failureShutdownDuration = failureShutdownDuration,
          writeBatchSize = writeBatchSize,
          writeThreshold = writeThreshold
        )
    }
  }
}
