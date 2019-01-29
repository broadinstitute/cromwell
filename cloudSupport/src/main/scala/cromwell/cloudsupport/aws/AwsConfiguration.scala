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

package cromwell.cloudsupport.aws

import cats.data.Validated._
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigException}
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.cloudsupport.aws.auth.{AssumeRoleMode, AwsAuthMode, CustomKeyMode, DefaultMode}
import net.ceedubs.ficus.Ficus._
import org.slf4j.LoggerFactory
import software.amazon.awssdk.regions.Region

final case class AwsConfiguration private (applicationName: String,
                                           authsByName: Map[String, AwsAuthMode],
                                           strRegion: Option[String]) {

  def auth(name: String): ErrorOr[AwsAuthMode] = {
    authsByName.get(name) match {
      case None =>
        val knownAuthNames = authsByName.keys.mkString(", ")
        s"`aws` configuration stanza does not contain an auth named '$name'.  Known auth names: $knownAuthNames".invalidNel
      case Some(a) => a.validNel
    }
  }

  def region: Option[Region] = strRegion.map(Region.of)
}

object AwsConfiguration {
  import scala.concurrent.duration._
  import scala.language.postfixOps

  lazy val DefaultConnectionTimeout = 3 minutes
  lazy val DefaultReadTimeout = 3 minutes

  private val log = LoggerFactory.getLogger("AwsConfiguration")

  final case class AwsConfigurationException(errorMessages: List[String]) extends MessageAggregation {
    override val exceptionContext = "AWS configuration"
  }

  def apply(config: Config): AwsConfiguration = {

    val awsConfig = config.getConfig("aws")

    val appName = validate { awsConfig.as[String]("application-name") }

    val region: Option[String] =
      awsConfig.getAs[String]("region")

    def buildAuth(authConfig: Config): ErrorOr[AwsAuthMode] = {

      def customKeyAuth(authConfig: Config, name: String, region: Option[String]): ErrorOr[AwsAuthMode] = validate {
        (authConfig.getAs[String]("access-key"), authConfig.getAs[String]("secret-key")) match {
          case (Some(accessKey), Some(secretKey)) =>
            CustomKeyMode(name, accessKey, secretKey, region)
          case _ => throw new ConfigException.Generic(s"""Access key and/or secret """ +
            s"""key missing for service account "$name". See reference.conf under the aws.auth, """ +
            s"""custom key section for details of required configuration.""")
        }
      }

      def defaultAuth(authConfig: Config, name: String, region: Option[String]): ErrorOr[AwsAuthMode] =  validate {
        DefaultMode(name, region)
      }

      def assumeRoleAuth(authConfig: Config, name: String, region: Option[String]): ErrorOr[AwsAuthMode] = validate {
        val externalId = authConfig.getOrElse("external-id", "")
        AssumeRoleMode(
          name,
          // We won't do anything with this now, but it is required for
          // assignment later
          authConfig.getString("base-auth"),
          authConfig.getString("role-arn"),
          externalId,
          region
        )
      }
      val name = authConfig.getString("name")
      val scheme = authConfig.getString("scheme")

      scheme match {
        case "default" => defaultAuth(authConfig, name, region)
        case "custom_keys" => customKeyAuth(authConfig, name, region)
        case "assume_role" => assumeRoleAuth(authConfig, name, region)
        case wut => s"Unsupported authentication scheme: $wut".invalidNel
      }
    }

    def assignDependency(dependentAuth: AssumeRoleMode, auths: List[ErrorOr[AwsAuthMode]]): Unit = {
      // We only care here about valid auth blocks. If something is invalid
      // we need to throw at some point anyway. This helps unwrap some of the
      // validation type wrappers that are involved at this point in the code
      val validAuths = auths.collect { case Valid(v) => v }

      // Look for the base auth from the config. If we find it, we'll assign
      // here. Unfortunately, we will rely on a runtime error if the base auth
      // does not end up getting assigned to the AssumeRoleMode object
      val baseAuth = validAuths.collectFirst { case a if a.name == dependentAuth.baseAuthName => a }
      baseAuth foreach dependentAuth.assign
    }

    def assignDependencies(auths: List[ErrorOr[AwsAuthMode]]): List[ErrorOr[AwsAuthMode]] = {
      // Assume role is somewhat special. We need to process assume role type
      // auths after its base auth is created. As such, we'll wire in the
      // base auth element after the list is created
      auths.collect { case Valid(arm: AssumeRoleMode) => assignDependency(arm, auths) }
      auths
    }
    val listOfErrorOrAuths: List[ErrorOr[AwsAuthMode]] =
      assignDependencies(awsConfig.as[List[Config]]("auths").map(buildAuth))
    val errorOrAuthList: ErrorOr[List[AwsAuthMode]] = listOfErrorOrAuths.sequence[ErrorOr, AwsAuthMode]

    def uniqueAuthNames(list: List[AwsAuthMode]): ErrorOr[Unit] = {
      val duplicateAuthNames = list.groupBy(_.name) collect { case (n, as) if as.size > 1 => n }
      if (duplicateAuthNames.nonEmpty) {
        ("Duplicate auth names: " + duplicateAuthNames.mkString(", ")).invalidNel
      } else {
        ().validNel
      }
    }

    (appName, errorOrAuthList).flatMapN { (name, list) =>
      uniqueAuthNames(list) map { _ =>
        AwsConfiguration(name, list map { a => a.name -> a } toMap, region)
      }
    } match {
      case Valid(r) => r
      case Invalid(f) =>
        val errorMessages = f.toList.mkString(", ")
        log.error(errorMessages)
        throw AwsConfigurationException(f.toList)
    }
  }
}
