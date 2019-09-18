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

import com.typesafe.config.{ConfigFactory}
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}

class AwsBatchConfigurationSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks with BeforeAndAfterAll {

  behavior of "AwsBatchConfigurationSpec"

  val mockFile = DefaultPathBuilder.createTempFile()

  override def afterAll(): Unit = {
    mockFile.delete(swallowIOExceptions = true)
    ()
  }

  val globalConfig = ConfigFactory.parseString(
    s"""
      |aws {
      |
      |  application-name = "cromwell"
      |   numSubmitAttempts = 6
      |   numCreateDefinitionAttempts = 6
      |
      |  auths = [
      |    {
      |      name = "application-default"
      |      scheme = "default"
      |    },
      |    {
      |      name = "user-via-refresh"
      |      scheme = "custom_keys"
      |      access-key = "secret_key"
      |      secret-key = "${mockFile.pathAsString}"
      |    },
      |    {
      |      name = "service-account"
      |      scheme = "default"
      |    }
      |  ]
      |}
      |
    """.stripMargin)

  val backendConfig = ConfigFactory.parseString(
    """
      |  // Base bucket for workflow executions
      |  root = "s3://my-cromwell-workflows-bucket"
      |
      |  auth = "application-default"
      |   numSubmitAttempts = 6
      |   numCreateDefinitionAttempts = 6
      |
      |  default-runtime-attributes {
      |      failOnStderr: false
      |      continueOnReturnCode: 0
      |      cpu: 1
      |      memory: "2 GB"
      |      disks: "local-disk"
      |      noAddress: false
      |      zones:["us-east-1a", "us-east-1b"]
      |  }
      |
      |  dockerhub {
      |    account = "dockerAccount"
      |    token = "dockerToken"
      |  }
      |
      |  filesystems {
      |    s3 {
      |      // A reference to a potentially different auth for manipulating files via engine functions.
      |      auth = "application-default"
      |    }
      |  }
      |
    """.stripMargin)

  it should "fail to instantiate if any required configuration is missing" in {

    val configs = Table(
      ("backendConfig", "globalConfig"),
      (backendConfig, globalConfig.withoutPath("aws")),
      (backendConfig.withoutPath("root"), globalConfig),
      (backendConfig.withoutPath("filesystems"), globalConfig),
      (backendConfig.withoutPath("filesystems.s3"), globalConfig),
      (backendConfig.withoutPath("filesystems.s3.auth"), globalConfig)
    )

    forAll(configs) { (backend, global) =>
      an[Exception] shouldBe thrownBy {
        new AwsBatchConfiguration(BackendConfigurationDescriptor(backend, global))
      }
    }
  }

  it should "have correct root" in {
    new AwsBatchConfiguration(BackendConfigurationDescriptor(backendConfig, globalConfig)).root shouldBe "s3://my-cromwell-workflows-bucket"
  }

  it should "have correct docker" in {
    val dockerConf = new AwsBatchConfiguration(BackendConfigurationDescriptor(backendConfig, globalConfig)).dockerCredentials
    dockerConf shouldBe defined
    dockerConf.get.token shouldBe "dockerToken"
  }
}
