package cromwell.engine.workflow.workflowstore

import java.util.UUID

import com.typesafe.config.{Config, ConfigFactory}
import io.circe.generic.auto._
import io.circe.literal._
import io.circe.syntax._
import io.circe._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps


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

  // An `Encoder[FiniteDuration]` is needed for the `.asJson.toString()` shenanigans below. Unfortunately the compiler
  // appears to get confused about this encoder actually being used. If this is made fully private or a local the compiler
  // wrongly emits a warning about it being unused, which our compiler flag settings then promote to an error.
  private [engine] implicit val finiteDurationEncoder: Encoder[FiniteDuration] = (d: FiniteDuration) => d.toString.asJson

  override def toString: String = this.asJson.toString()
}

object WorkflowHeartbeatConfig {

  def apply(config: Config): WorkflowHeartbeatConfig = {
    val cromwellId: String = config.as[Option[String]]("system.cromwell_id").getOrElse("cromid-" + UUID.randomUUID().toString.take(7))

    val heartbeats: Config = config.getOrElse("system.workflow-heartbeats", ConfigFactory.empty())
    // Default to 10 minutes and don't allow values less than 10 seconds even for testing purposes.
    val minHeartbeatTtl = 10 seconds
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
