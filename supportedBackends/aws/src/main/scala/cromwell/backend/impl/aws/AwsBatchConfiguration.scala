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

import cromwell.filesystems.s3.S3PathBuilderFactory
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.core.{BackendDockerConfiguration}
import cromwell.core.path.PathBuilderFactory
import cromwell.cloudsupport.aws.AwsConfiguration

class AwsBatchConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val awsConfig = AwsConfiguration(configurationDescriptor.globalConfig)

  val root = configurationDescriptor.backendConfig.getString("root")
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig
  val batchAttributes = AwsBatchAttributes.fromConfigs(awsConfig, configurationDescriptor.backendConfig)
  val awsAuth = batchAttributes.auth
  val dockerCredentials = BackendDockerConfiguration.build(configurationDescriptor.backendConfig).dockerCredentials
  val fileSystem =
    configurationDescriptor.backendConfig.hasPath("filesystems.s3") match {
      case true =>  "s3"
      case false => "local"
  }
  val pathBuilderFactory = configurationDescriptor.backendConfig.hasPath("filesystems.s3") match {
    case true => S3PathBuilderFactory(configurationDescriptor.globalConfig, configurationDescriptor.backendConfig)
    case false =>
    PathBuilderFactory
  }
}

object AWSBatchStorageSystems {
  val s3:String = "s3"
  val efs:String = "efs"
  val ebs:String = "ebs"
  val local:String = "local"
}

