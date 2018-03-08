package cromwell.core

import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration.FiniteDuration

object LoadConfig {
  private val conf = ConfigFactory.load().getConfig("load-control")
  val JobStoreReadThreshold = conf.as[Int]("job-store-read")
  val JobStoreWriteThreshold = conf.as[Int]("job-store-write")
  val CallCacheReadThreshold = conf.as[Int]("call-cache-read")
  val CallCacheWriteThreshold = conf.as[Int]("call-cache-write")
  val KeyValueReadThreshold = conf.as[Int]("key-value-read")
  val KeyValueWriteThreshold = conf.as[Int]("key-value-write")
  val IoThreshold = conf.as[Int]("io")
  val MetadataWriteThreshold = conf.as[Int]("metadata-write")
  val MonitoringFrequency = conf.as[FiniteDuration]("monitoring-frequency")
  val MemoryThresholdInMB = conf.as[Int]("memory-threshold-in-mb")
  val MemoryMeasurementWindow = conf.as[Int]("memory-measurement-window")
  val PAPIThreshold = conf.as[Int]("papi-requests")
}
