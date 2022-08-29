/*
 * Copyright 2020 Amazon.com, Inc. or its affiliates.
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

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.typesafe.config.Config
import cromwell.services.metadata.MetadataService.{PutMetadataAction, PutMetadataActionAndRespond}
import cromwell.services.metadata.impl.MetadataServiceActor


/**
  * A metadata service implementation which will function as a standard metadata service but also push all metadata
  * events to AWS SNS (Simple Notification Service). This class closely follows the pattern established in the
  * HybridPubSubMetadataServiceActor
  *
  * Under the hood it maintains its own MetadataServiceActor and AwsSnsMetadataServiceActor. All messages are routed
  * to the MetadataServiceActor. PutMetadataActions are also sent to the AwsSnsMetadataServiceActor. PutMetadataActionAndRespond
  * messages will be sent to the SnsMetadataServiceActor as a standard PutMetadataAction, i.e. only the standard
  * metadata service will be ACKing the request.
  *
  * To use this actor something similar to the following should be present in the cromwell.conf file:
  * <pre>
  * services {
  *   MetadataService {
  *     class="cromwell.services.metadata.impl.sns.HybridSnsMetadataServiceActor"
  *     config {
  *       aws {
  *         application-name = "cromwell"
  *         auths = [{
  *           name = "default"
  *           scheme = "default"
  *         }]
  *         region = "us-east-1"
  *         topicArn = "arn:aws:sns:us-east-1:1111111111111:cromwell-metadata"
  *       }
  *     }
  *   }
  * }
  * </pre>
  *
  * @see cromwell.services.metadata.impl.sns.AwsSnsMetadataServiceActor
  */
class HybridSnsMetadataServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with ActorLogging  {
  val standardMetadataActor: ActorRef = context.actorOf(MetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))
  val awsSnsMetadataActor: ActorRef = context.actorOf(AwsSnsMetadataServiceActor.props(serviceConfig, globalConfig, serviceRegistryActor))

  override def receive = {
    case action: PutMetadataAction =>
      standardMetadataActor forward action
      awsSnsMetadataActor forward action
    case action: PutMetadataActionAndRespond =>
      standardMetadataActor forward action
      awsSnsMetadataActor forward PutMetadataAction(action.events)
    case anythingElse => standardMetadataActor forward anythingElse
  }
}

