package cromwell.core

import org.scalatest.{FlatSpec, Matchers}
import scala.concurrent.duration._

class LoadConfigSpec extends FlatSpec with Matchers {
  behavior of "LoadConfig"
  
  it should "parse load config" in {
    LoadConfig.JobStoreReadThreshold shouldBe 10000
    LoadConfig.JobStoreWriteThreshold shouldBe 10000
    LoadConfig.CallCacheReadThreshold shouldBe 1000
    LoadConfig.CallCacheWriteThreshold shouldBe 10000
    LoadConfig.KeyValueReadThreshold shouldBe 10000
    LoadConfig.KeyValueWriteThreshold shouldBe 10000
    LoadConfig.MetadataWriteThreshold shouldBe 100000
    LoadConfig.IoThreshold shouldBe 10000
    LoadConfig.MonitoringFrequency shouldBe 5.seconds
    LoadConfig.MemoryThresholdInMB shouldBe 1024
    LoadConfig.PAPIThreshold shouldBe 10000
    LoadConfig.MemoryMeasurementWindow shouldBe 6
  }
}
