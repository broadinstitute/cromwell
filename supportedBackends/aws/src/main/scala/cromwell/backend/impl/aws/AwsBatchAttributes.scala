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

import java.net.{URI, URL}

import cats.data.Validated._
import com.typesafe.config.{Config, ConfigValue}
import cromwell.cloudsupport.aws.auth.AwsAuthMode
import cromwell.backend.impl.aws.callcaching.{CopyCachedOutputs, AwsBatchCacheHitDuplicationStrategy, UseOriginalCachedOutputs}
import cromwell.cloudsupport.aws.AwsConfiguration
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.{StringReader, ValueReader}
import org.slf4j.{Logger, LoggerFactory}

import scala.collection.JavaConverters._

case class AwsBatchAttributes(auth: AwsAuthMode,
                         executionBucket: String,
                         duplicationStrategy: AwsBatchCacheHitDuplicationStrategy)

object AwsBatchAttributes {
  lazy val Logger = LoggerFactory.getLogger("AwsBatchAttributes")

  private val availableConfigKeys = Set(
    "concurrent-job-limit",
    "root",
    "dockerhub",
    "dockerhub.account",
    "dockerhub.token",
    "filesystems",
    "filesystems.s3.auth",
    "filesystems.s3.caching.duplication-strategy",
    "default-runtime-attributes",
    "default-runtime-attributes.disks",
    "default-runtime-attributes.memory",
    "default-runtime-attributes.zones",
    "default-runtime-attributes.continueOnReturnCode",
    "default-runtime-attributes.cpu",
    "default-runtime-attributes.noAddress",
    "default-runtime-attributes.docker",
    "default-runtime-attributes.queueArn",
    "default-runtime-attributes.failOnStderr"
  )

  private val deprecatedAwsBatchKeys: Map[String, String] = Map(
  )

  private val context = "AwsBatch"

  implicit val urlReader: ValueReader[URL] = StringReader.stringValueReader.map { URI.create(_).toURL }

  def apply(awsConfig: AwsConfiguration, backendConfig: Config): AwsBatchAttributes = {
    val configKeys = backendConfig.entrySet().asScala.toSet map { entry: java.util.Map.Entry[String, ConfigValue] => entry.getKey }
    warnNotRecognized(configKeys, availableConfigKeys, context, Logger)

    def warnDeprecated(keys: Set[String], deprecated: Map[String, String], context: String, logger: Logger) = {
      val deprecatedKeys = keys.intersect(deprecated.keySet)
      deprecatedKeys foreach { key => logger.warn(s"Found deprecated configuration key $key, replaced with ${deprecated.get(key)}") }
    }

    warnDeprecated(configKeys, deprecatedAwsBatchKeys, context, Logger)

    val executionBucket: ErrorOr[String] = validate { backendConfig.as[String]("root") }
    val filesystemAuthName: ErrorOr[String] = validate { backendConfig.as[String]("filesystems.s3.auth") }
    val duplicationStrategy = validate { backendConfig.as[Option[String]]("filesystems.s3.caching.duplication-strategy").getOrElse("copy") match {
      case "copy" => CopyCachedOutputs
      case "reference" => UseOriginalCachedOutputs
      case other => throw new IllegalArgumentException(s"Unrecognized caching duplication strategy: $other. Supported strategies are copy and reference. See reference.conf for more details.")
    } }


    def authAwsConfigForAwsBatchAttributes(
                          bucket: String,
                          authName: String,
                          cachingStrategy: AwsBatchCacheHitDuplicationStrategy): ErrorOr[AwsBatchAttributes] = validate {
      awsConfig.auth(authName) match {
        case Valid(a) => AwsBatchAttributes(a, bucket, cachingStrategy)
        case Invalid(a) => throw new IllegalArgumentException(s"AwsBatch Configuration, auth name '$a' is not valid.")
      }
    }

    (executionBucket, filesystemAuthName, duplicationStrategy) flatMapN authAwsConfigForAwsBatchAttributes match {
      case Valid(r) => r
      case Invalid(f) =>
        throw new IllegalArgumentException with MessageAggregation {
          override val exceptionContext = "AwsBatch Configuration is not valid: Errors"
          override val errorMessages = f.toList
        }
    }
  }
}
