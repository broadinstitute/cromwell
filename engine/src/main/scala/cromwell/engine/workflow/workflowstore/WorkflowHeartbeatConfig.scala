package cromwell.engine.workflow.workflowstore

import java.util.UUID

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.engine.workflow.workflowstore.WorkflowHeartbeatConfig._
import io.circe._
import io.circe.generic.semiauto._
import io.circe.syntax._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps
import mouse.all._

/**
  *
  * @param cromwellId         Identifier for this Cromwell for horizontaling purposes.
  * @param heartbeatInterval  How long between `WorkflowActor` heartbeats and `WorkflowStoreHeartbeatWriteActor` database flushes.
  * @param ttl                Entries in the workflow store with heartbeats older than `ttl` ago are presumed abandoned and up for grabs.
  * @param writeBatchSize     The maximum size of a write batch.
  * @param writeThreshold     The threshold of heartbeat writes above which load is considered high.
  */
case class WorkflowHeartbeatConfig(
                                    cromwellId: String,
                                    heartbeatInterval: FiniteDuration,
                                    ttl: FiniteDuration,
                                    writeBatchSize: Int,
                                    writeThreshold: Int)
{
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
    val randomSuffix = config
      .getOrElse("system.cromwell_id_random_suffix", true)
      .fold("-" + UUID.randomUUID().toString.take(7), "")
    val cromwellId: String = config.getOrElse("system.cromwell_id", "cromid") + randomSuffix

    val heartbeats: Config = config.getOrElse("system.workflow-heartbeats", ConfigFactory.empty())
    // Default to 10 mins and don't allow values less than 10 secs except for testing purposes.
    val minHeartbeatTtl = config.getOrElse("danger.debug.only.minimum-heartbeat-ttl", 10.seconds)
    val ttl: FiniteDuration = heartbeats.getOrElse("ttl", 10 minutes).max(minHeartbeatTtl)
    // Default to one third the TTL and don't allow values less than one third of the minimum heartbeat TTL.
    val heartbeatInterval: FiniteDuration = heartbeats.getOrElse("heartbeat-interval", ttl / 3).max(minHeartbeatTtl / 3)

    // Default and sanity bound all other parameters as well.
    val writeBatchSize: Int = Math.max(heartbeats.getOrElse("write-batch-size", 100), 1)
    val writeThreshold: Int = Math.max(heartbeats.getOrElse("write-threshold", 100), 1)

    WorkflowHeartbeatConfig(
      cromwellId = cromwellId,
      heartbeatInterval = heartbeatInterval,
      ttl = ttl,
      writeBatchSize = writeBatchSize,
      writeThreshold = writeThreshold
    )
  }
}
