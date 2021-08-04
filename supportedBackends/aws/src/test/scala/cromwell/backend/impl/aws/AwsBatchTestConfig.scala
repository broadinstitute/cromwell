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

import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendConfigurationDescriptor

object AwsBatchTestConfig {

  private val AwsBatchBackendConfigString =
    """
      |root = "s3://my-cromwell-workflows-bucket"
      |
      |filesystems {
      |  s3 {
      |    auth = "default"
      |  }
      |}
      |
      |auth = "default"
      |   numSubmitAttempts = 6
      |   numCreateDefinitionAttempts = 6
      |
      |default-runtime-attributes {
      |    cpu: 1
      |    failOnStderr: false
      |    continueOnReturnCode: 0
      |    docker: "ubuntu:latest"
      |    memory: "2 GB"
      |    disks: "local-disk"
      |    noAddress: false
      |    zones:["us-east-1a", "us-east-1b"]
      |    queueArn: "arn:aws:batch:us-east-1:111222333444:job-queue/job-queue"
      |    scriptBucketName: "my-bucket"
      |    awsBatchRetryAttempts: 1
      |}
      |
      |""".stripMargin

  private val NoDefaultsConfigString =
    """
      |root = "s3://my-cromwell-workflows-bucket"
      |
      |auth = "default"
      |   numSubmitAttempts = 6
      |   numCreateDefinitionAttempts = 6
      |
      |filesystems {
      |  s3 {
      |    auth = "default"
      |  }
      |}
      |""".stripMargin

  private val AwsBatchGlobalConfigString =
    s"""
       |aws {
       |  application-name = "cromwell"
       |  auths = [
       |    {
       |      name = "default"
       |      scheme = "default"
       |    }
       |  ]
       |}
       |
       |backend {
       |  default = "AWS"
       |  providers {
       |    AWS {
       |      actor-factory = "cromwell.backend.impl.aws.AwsBatchBackendLifecycleFactory"
       |      config {
       |      $AwsBatchBackendConfigString
       |      }
       |    }
       |  }
       |}
       |
       |""".stripMargin

  val AwsBatchBackendConfig = ConfigFactory.parseString(AwsBatchBackendConfigString)
  val AwsBatchGlobalConfig = ConfigFactory.parseString(AwsBatchGlobalConfigString)
  val AwsBatchBackendNoDefaultConfig = ConfigFactory.parseString(NoDefaultsConfigString)
  val AwsBatchBackendConfigurationDescriptor = BackendConfigurationDescriptor(AwsBatchBackendConfig, AwsBatchGlobalConfig)
  val NoDefaultsConfigurationDescriptor = BackendConfigurationDescriptor(AwsBatchBackendNoDefaultConfig, AwsBatchGlobalConfig)
}

object AwsBatchTestConfigForLocalFS {

  private val AwsBatchBackendConfigString =
    """
      |root = "/root/my-cromwell-workflows"
      |
      |filesystems {
      |  local {
      |    auth = "default"
      |  }
      |}
      |
      |auth = "default"
      |   numSubmitAttempts = 6
      |   numCreateDefinitionAttempts = 6
      |
      |default-runtime-attributes {
      |    cpu: 1
      |    failOnStderr: false
      |    continueOnReturnCode: 0
      |    docker: "ubuntu:latest"
      |    memory: "2 GB"
      |    disks: "local-disk"
      |    noAddress: false
      |    zones:["us-east-1a", "us-east-1b"]
      |    queueArn: "arn:aws:batch:us-east-1:111222333444:job-queue/job-queue"
      |    scriptBucketName: ""
      |    awsBatchRetryAttempts: 1
      |}
      |
      |""".stripMargin

  private val NoDefaultsConfigString =
    """
      |root = "/root/my-cromwell-workflows"
      |
      |auth = "default"
      |   numSubmitAttempts = 6
      |   numCreateDefinitionAttempts = 6
      |
      |filesystems {
      |  local {
      |    auth = "default"
      |  }
      |}
      |""".stripMargin

  private val AwsBatchGlobalConfigString =
    s"""
       |aws {
       |  application-name = "cromwell"
       |  auths = [
       |    {
       |      name = "default"
       |      scheme = "default"
       |    }
       |  ]
       |}
       |
       |backend {
       |  default = "AWS"
       |  providers {
       |    AWS {
       |      actor-factory = "cromwell.backend.impl.aws.AwsBatchBackendLifecycleFactory"
       |      config {
       |      $AwsBatchBackendConfigString
       |      }
       |    }
       |  }
       |}
       |
       |""".stripMargin

  val AwsBatchBackendConfig = ConfigFactory.parseString(AwsBatchBackendConfigString)
  val AwsBatchGlobalConfig = ConfigFactory.parseString(AwsBatchGlobalConfigString)
  val AwsBatchBackendNoDefaultConfig = ConfigFactory.parseString(NoDefaultsConfigString)
  val AwsBatchBackendConfigurationDescriptor = BackendConfigurationDescriptor(AwsBatchBackendConfig, AwsBatchGlobalConfig)
  val NoDefaultsConfigurationDescriptor = BackendConfigurationDescriptor(AwsBatchBackendNoDefaultConfig, AwsBatchGlobalConfig)
}