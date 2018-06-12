/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
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

package cromwell.backend.impl.aws

import java.util.UUID

import akka.actor.Props
import akka.testkit._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendWorkflowInitializationActor.{InitializationFailed, Initialize}
import cromwell.backend.async.RuntimeAttributeValidationFailures
import cromwell.backend.{BackendConfigurationDescriptor, BackendSpec, BackendWorkflowDescriptor}
// import cromwell.cloudsupport.aws.auth.AwsAuthModeSpec
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.Tags.PostWomTest
import cromwell.core.{TestKitSuite}
// import cromwell.core.logging.LoggingTest._
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import spray.json._
import wom.graph.CommandCallNode

import scala.concurrent.duration._

class AwsBatchInitializationActorSpec extends TestKitSuite("AwsBatchInitializationActorSpec") with FlatSpecLike with Matchers
  with ImplicitSender with Mockito {
  val Timeout: FiniteDuration = 10.second.dilated

  import BackendSpec._

  val HelloWorld: String =
    s"""
      |task hello {
      |  String addressee = "you"
      |  command {
      |    echo "Hello $${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow wf_hello {
      |  call hello
      |}
    """.stripMargin

  val globalConfig: Config = ConfigFactory.parseString(
    """
      |aws {
      |
      |  application-name = "cromwell"
      |
      |  auths = [
      |    {
      |      name = "application-default"
      |      scheme = "application_default"
      |    },
      |    {
      |      name = "user-via-refresh"
      |      scheme = "refresh_token"
      |      access-key = "secret_id"
      |      secret-key = "secret_secret"
      |    }
      |  ]
      |}
      |""".stripMargin)

  val backendConfigTemplate: String =
    """
      |  // Base bucket for workflow executions
      |  root = "s3://my-cromwell-workflows-bucket"
      |
      |  // Polling for completion backs-off gradually for slower-running jobs.
      |  // This is the maximum polling interval (in seconds):
      |  maximum-polling-interval = 600
      |
      |  default-runtime-attributes {
      |     cpu: 1
      |     failOnStderr: false
      |     # Allowed to be a boolean, or a list of Ints, or an Int
      |     continueOnReturnCode: 0
      |     memory: "2 GB"
      |     # Allowed to be a String, or a list of Strings
      |     disks: "local-disk"
      |     noAddress: false
      |     zones: ["us-east-1a", "us-east-1b"]
      |  }
      |  filesystems {
      |    s3 {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "application-default"
      |    }
      |  }
      |
      |[DOCKERHUBCONFIG]
      |""".stripMargin

  val refreshTokenConfigTemplate: String =
    """
      |  // Base bucket for workflow executions
      |  root = "s3://my-cromwell-workflows-bucket"
      |
      |  // Polling for completion backs-off gradually for slower-running jobs.
      |  // This is the maximum polling interval (in seconds):
      |  maximum-polling-interval = 600
      |
      |  default-runtime-attributes {
      |     cpu: 1
      |     failOnStderr: false
      |     # Allowed to be a boolean, or a list of Ints, or an Int
      |     continueOnReturnCode: 0
      |     memory: "2 GB"
      |     # Allowed to be a String, or a list of Strings
      |     disks: "local-disk"
      |     noAddress: false
      |     zones: ["us-east-1a", "us-east-1b"]
      |  }
      |
      |  filesystems {
      |    s3 {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "user-via-refresh"
      |    }
      |  }
      |""".stripMargin

  val backendConfig: Config = ConfigFactory.parseString(backendConfigTemplate.replace("[DOCKERHUBCONFIG]", ""))

  val dockerBackendConfig: Config = ConfigFactory.parseString(backendConfigTemplate.replace("[DOCKERHUBCONFIG]",
    """
      |dockerhub {
      |  account = "my@docker.account"
      |  token = "mydockertoken"
      |}
      | """.stripMargin))

  val defaultBackendConfig = BackendConfigurationDescriptor(backendConfig, globalConfig)

  val refreshTokenConfig: Config = ConfigFactory.parseString(refreshTokenConfigTemplate)

  private def getAwsBatchBackendProps(workflowDescriptor: BackendWorkflowDescriptor,
                                 calls: Set[CommandCallNode],
                                 configuration: AwsBatchConfiguration): Props = {
    val ioActor = mockIoActor
    val params = AwsBatchInitializationActorParams(workflowDescriptor, ioActor, calls, configuration, emptyActor, restarting = false)
    Props(new AwsBatchInitializationActor(params)).withDispatcher(BackendDispatcher)
  }

  private def getAwsBatchBackend(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[CommandCallNode], conf: BackendConfigurationDescriptor) = {
    val props = getAwsBatchBackendProps(workflowDescriptor, calls, new AwsBatchConfiguration(conf))
    system.actorOf(props, "TestableAwsBatchInitializationActor-" + UUID.randomUUID)
  }

  behavior of "AwsBatchInitializationActor"

  // it should "log a warning message when there are unsupported runtime attributes" taggedAs IntegrationTest in {
  //   AwsAuthModeSpec.assumeHasApplicationDefaultCredentials()
  //
  //   within(Timeout) {
  //     val workflowDescriptor = buildWdlWorkflowDescriptor(HelloWorld,
  //       runtime = """runtime { docker: "ubuntu/latest" test: true }""")
  //     val backend = getAwsBatchBackend(workflowDescriptor, workflowDescriptor.callable.taskCallNodes,
  //       defaultBackendConfig)
  //     val eventPattern =
  //       "Key/s [test] is/are not supported by backend. Unsupported attributes will not be part of job executions."
  //     EventFilter.warning(pattern = escapePattern(eventPattern), occurrences = 1) intercept {
  //       backend ! Initialize
  //     }
  //     expectMsgPF() {
  //       case InitializationSuccess(_) => //Docker entry is present.
  //       case InitializationFailed(failure) => fail(s"InitializationSuccess was expected but got $failure")
  //     }
  //   }
  // }

  // Depends on https://github.com/broadinstitute/cromwell/issues/2606
  it should "return InitializationFailed when docker runtime attribute key is not present" taggedAs PostWomTest ignore {
    within(Timeout) {
      val workflowDescriptor = buildWdlWorkflowDescriptor(HelloWorld, runtime = """runtime { }""")
      val backend = getAwsBatchBackend(workflowDescriptor, workflowDescriptor.callable.taskCallNodes,
        defaultBackendConfig)
      backend ! Initialize
      expectMsgPF() {
        case InitializationFailed(failure) =>
          failure match {
            case exception: RuntimeAttributeValidationFailures =>
              if (!exception.getMessage.equals("Runtime validation failed:\nTask hello has an invalid runtime attribute docker = !! NOT FOUND !!"))
                fail("Exception message is not equal to 'Runtime validation failed:\nTask hello has an invalid runtime attribute docker = !! NOT FOUND !!'.")
          }
      }
    }
  }

  private case class TestingBits(actorRef: TestActorRef[AwsBatchInitializationActor], configuration: AwsBatchConfiguration)

}

object AwsBatchInitializationActorSpec {
  def normalize(str: String) = {
    str.parseJson.prettyPrint
  }
}
