package cromwell.engine.workflow.workflowstore

import com.typesafe.config.ConfigFactory
import common.exception.AggregatedMessageException
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._

class WorkflowHeartbeatConfigSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "WorkflowHeartbeatConfig"

  it should "generate a default config" in {
    val workflowHeartbeatConfig = WorkflowHeartbeatConfig(WorkflowHeartbeatConfigSpec.DefaultConfig)
    workflowHeartbeatConfig.cromwellId should startWith("cromid-")
    workflowHeartbeatConfig.cromwellId should have length 14
    workflowHeartbeatConfig.heartbeatInterval should be (2.minutes)
    workflowHeartbeatConfig.ttl should be(10.minutes)
    workflowHeartbeatConfig.failureShutdownDuration should be(5.minutes)
    workflowHeartbeatConfig.writeBatchSize should be(10000)
    workflowHeartbeatConfig.writeThreshold should be(10000)
  }

  private val validConfigTests = Table(
    ("description", "config", "expected"),
    (
      "from an empty config",
      "system.cromwell_id_random_suffix = false",
      WorkflowHeartbeatConfig("cromid", 2.minutes, 10.minutes, 5.minutes, 10000, 10000),
    ),
    (
      "with a specified cromid",
      """system.cromwell_id = "new_crom_name"""",
      WorkflowHeartbeatConfig("new_crom_name", 2.minutes, 10.minutes, 5.minutes, 10000, 10000),
    ),
    (
      "with a specified heartbeat interval",
      "system.workflow-heartbeats.heartbeat-interval = 3 minutes",
      WorkflowHeartbeatConfig("cromid", 3.minutes, 10.minutes, 5.minutes, 10000, 10000),
    ),
    (
      "with a specified ttl",
      "system.workflow-heartbeats.ttl = 5 minutes",
      WorkflowHeartbeatConfig("cromid", 2.minutes, 5.minutes, 5.minutes, 10000, 10000),
    ),
    (
      "with a ttl less than the default heartbeat interval",
      """|system.workflow-heartbeats.ttl = 1 minutes
         |system.workflow-heartbeats.heartbeat-interval = 59 seconds
         |system.workflow-heartbeats.write-failure-shutdown-duration = 0 minutes
         |""".stripMargin,
      WorkflowHeartbeatConfig("cromid", 59.seconds, 1.minutes, 0.minutes, 10000, 10000),
    ),
    (
      "with a specified shutdown duration",
      "system.workflow-heartbeats.write-failure-shutdown-duration = 1 minute",
      WorkflowHeartbeatConfig("cromid", 2.minutes, 10.minutes, 1.minute, 10000, 10000),
    ),
    (
      "with a specified batch size",
      "system.workflow-heartbeats.write-batch-size = 2000",
      WorkflowHeartbeatConfig("cromid", 2.minutes, 10.minutes, 5.minutes, 2000, 10000),
    ),
    (
      "with a specified threshold",
      "system.workflow-heartbeats.write-threshold = 5000",
      WorkflowHeartbeatConfig("cromid", 2.minutes, 10.minutes, 5.minutes, 10000, 5000),
    ),
    (
      "when trying to set the ttl below the minimum",
      """|system.workflow-heartbeats.ttl = 9 seconds
         |system.workflow-heartbeats.heartbeat-interval = 8 seconds
         |system.workflow-heartbeats.write-failure-shutdown-duration = 0 minutes
         |""".stripMargin,
      WorkflowHeartbeatConfig("cromid", 8.seconds, 10.seconds, 0.minutes, 10000, 10000),
    ),
    (
      "when trying to set the interval below the minimum",
      "system.workflow-heartbeats.heartbeat-interval = 3 seconds",
      WorkflowHeartbeatConfig("cromid", 10.seconds / 3, 10.minutes, 5.minutes, 10000, 10000),
    ),
    (
      "when trying to set a negative shutdown duration",
      "system.workflow-heartbeats.write-failure-shutdown-duration = -1 seconds",
      WorkflowHeartbeatConfig("cromid", 2.minutes, 10.minutes, 0.minutes, 10000, 10000),
    ),
  )

  forAll(validConfigTests) { (description, configString, expected) =>
    it should s"create an instance $description" in {
      val config = ConfigFactory.parseString(
        // Remove the randomness from the cromid
        s"""|system.cromwell_id_random_suffix = false
            |$configString
            |""".stripMargin
      ).withFallback(WorkflowHeartbeatConfigSpec.DefaultConfig)
      WorkflowHeartbeatConfig(config) should be(expected)
    }
  }

  private val invalidConfigTests = Table(
    ("description", "config", "expected"),
    (
      "when the heartbeat interval is equal to the ttl",
      """|system.workflow-heartbeats {
         |  ttl = 2 minutes
         |  write-failure-shutdown-duration = 0 seconds
         |}
         |""".stripMargin,
      AggregatedMessageException(
        "Errors parsing WorkflowHeartbeatConfig",
        List(
          "The system.workflow-heartbeats.heartbeat-interval (2 minutes)" +
            " is not less than the system.workflow-heartbeats.ttl (2 minutes).",
        )
      )
    ),
    (
      "when the heartbeat interval is not less than the ttl",
      """|system.workflow-heartbeats {
         |  ttl = 1 minute
         |  write-failure-shutdown-duration = 0 seconds
         |}
         |""".stripMargin,
      AggregatedMessageException(
        "Errors parsing WorkflowHeartbeatConfig",
        List(
          "The system.workflow-heartbeats.heartbeat-interval (2 minutes)" +
            " is not less than the system.workflow-heartbeats.ttl (1 minute).",
        )
      )
    ),
    (
      "when the failure duration is greater than the ttl",
      """|system.workflow-heartbeats {
         |  ttl = 5 minutes
         |  write-failure-shutdown-duration = 301 seconds
         |}
         |""".stripMargin,
      AggregatedMessageException(
        "Errors parsing WorkflowHeartbeatConfig",
        List(
          "The system.workflow-heartbeats.write-failure-shutdown-duration (301 seconds)" +
            " is greater than the system.workflow-heartbeats.ttl (5 minutes).",
        )
      )
    ),
  )

  forAll(invalidConfigTests) { (description, configString, expected: AggregatedMessageException) =>
    it should s"fail to create an instance $description" in {
      val config = ConfigFactory.parseString(configString).withFallback(WorkflowHeartbeatConfigSpec.DefaultConfig)
      val exception = the[AggregatedMessageException] thrownBy WorkflowHeartbeatConfig(config)
      exception.getMessage should be(expected.getMessage)
      exception.errorMessages should contain theSameElementsAs expected.errorMessages
    }
  }
}

object WorkflowHeartbeatConfigSpec {
  val DefaultConfig = ConfigFactory.load()
}
