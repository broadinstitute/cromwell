package cromwell.services.metadata.impl.aws

import java.time.OffsetDateTime

import akka.actor.{ActorInitializationException, ActorRef, Props}
import akka.testkit.{EventFilter, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
import cromwell.services.ServicesSpec
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}


class AwsSnsMetadataServiceActorSpec extends ServicesSpec {
  import AwsSnsMetadataServiceActorSpec._

  val registryProbe: ActorRef = TestProbe().ref

  "An AwsSnsMetadataActor with an empty serviceConfig" should {
    "fail to build" in {
      EventFilter[ActorInitializationException](occurrences = 1) intercept {
        system.actorOf(Props(new AwsSnsMetadataServiceActor(emptyConfig, emptyConfig, registryProbe)))
      }
    }
  }

  "An AwsSnsMetadataActor with a topic and configuration" should {
    "successfully build" in {
      system.actorOf(Props(new AwsSnsMetadataServiceActor(configWithTopic, emptyConfig, registryProbe)))
    }

    "process an event" in {
      val actor = system.actorOf(Props(new AwsSnsMetadataServiceActor(configWithTopic, emptyConfig, registryProbe)))
      actor ! event
    }
  }
}

object AwsSnsMetadataServiceActorSpec {

  // This doesn't include a topic so should be a failure
  val emptyConfig: Config = ConfigFactory.empty()

  val configWithTopic: Config = ConfigFactory.parseString(
    """
      |aws {
      |  application-name = "cromwell"
      |  auths = [{
      |    name = "default"
      |    scheme = "default"
      |  }]
      |  region = "us-east-1"
      |  topicArn = "arn:aws:sns:us-east-1:1111111111111:cromwell-metadata"
      |}
    """.stripMargin
  )

  val event: MetadataEvent = MetadataEvent(MetadataKey(WorkflowId.randomId(), None, "key"),
    Option(MetadataValue("value")), OffsetDateTime.now)
}

