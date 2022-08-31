/*
 * Copyright 2022 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

package cromwell.services.metadata.impl.aws

import java.time.OffsetDateTime

import akka.actor.{ActorInitializationException, ActorRef, Props}
import akka.testkit.{EventFilter, TestProbe}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
import cromwell.services.ServicesSpec
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}


class AwsEventBridgeMetadataServiceActorSpec extends ServicesSpec {
  import AwsEventBridgeMetadataServiceActorSpec._

  val registryProbe: ActorRef = TestProbe().ref

  "An AwsEventBridgeMetadataActor with an empty serviceConfig" should {
    "fail to build" in {
      EventFilter[ActorInitializationException](occurrences = 1) intercept {
        system.actorOf(Props(new AwsEventBridgeMetadataServiceActor(emptyConfig, emptyConfig, registryProbe)))
      }
    }
  }

  "An AwsEventBridgeMetadataActor with a bus name and configuration" should {
    "successfully build" in {
      system.actorOf(Props(new AwsEventBridgeMetadataServiceActor(configWithBus, emptyConfig, registryProbe)))
    }

    "process an event" in {
      val actor = system.actorOf(Props(new AwsEventBridgeMetadataServiceActor(configWithBus, emptyConfig, registryProbe)))
      actor ! event
    }
  }
}

object AwsEventBridgeMetadataServiceActorSpec {

  // This doesn't include a topic so should be a failure
  val emptyConfig: Config = ConfigFactory.empty()

  val configWithBus: Config = ConfigFactory.parseString(
    """
      |aws {
      |  application-name = "cromwell"
      |  auths = [{
      |    name = "default"
      |    scheme = "default"
      |  }]
      |  region = "us-east-1"
      |  busName = "cromwell-metadata"
      |}
    """.stripMargin
  )

  val event: MetadataEvent = MetadataEvent(MetadataKey(WorkflowId.randomId(), None, "key"),
    Option(MetadataValue("value")), OffsetDateTime.now)
}

