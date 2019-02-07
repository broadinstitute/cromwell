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
package cromwell.cloudsupport.aws.auth

import com.google.api.client.json.jackson2.JacksonFactory
import cromwell.cloudsupport.aws.auth.AwsAuthMode.OptionLookup
import org.slf4j.LoggerFactory
import software.amazon.awssdk.auth.credentials._
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.sts.StsClient
import software.amazon.awssdk.services.sts.model.{AssumeRoleRequest, GetCallerIdentityRequest}

import scala.util.{Failure, Success, Try}

object AwsAuthMode {
  type OptionLookup = String => String
  lazy val jsonFactory = JacksonFactory.getDefaultInstance
}

sealed trait AwsAuthMode {
  protected lazy val log = LoggerFactory.getLogger(getClass.getSimpleName)

  /**
    * Validate the auth mode against provided options
    */
  def validate(options: OptionLookup): Unit = {
    ()
  }

  def name: String

  def credential(options: OptionLookup): AwsCredentials

  /**
    * Enables swapping out credential validation for various testing purposes ONLY.
    *
    * All traits in this file are sealed, all classes final, meaning things
    * like Mockito or other java/scala overrides cannot work.
    */
   private[auth] var credentialValidation: (AwsCredentials, Option[String]) => Unit =
     (credentials: AwsCredentials, region: Option[String]) => {
       val builder = StsClient.builder

       //If the region argument exists in config, set it in the builder.
       //Otherwise it is left unset and the AwsCredential builder will look in various places to supply,
       //ultimately using US-EAST-1 if none is found
       region.map(Region.of).foreach(builder.region)

       builder.credentialsProvider(StaticCredentialsProvider.create(credentials))
         .build
         .getCallerIdentity(GetCallerIdentityRequest.builder.build)
       ()
     }

  protected def validateCredential(credential: AwsCredentials, region: Option[String]) = {
    Try(credentialValidation(credential, region)) match {
      case Failure(ex) => throw new RuntimeException(s"Credentials are invalid: ${ex.getMessage}", ex)
      case Success(_) => credential
    }
  }
}

case object MockAuthMode extends AwsAuthMode {
  override val name = "no_auth"

  lazy val _credential = AnonymousCredentialsProvider.create.resolveCredentials()

  override def credential(options: OptionLookup): AwsCredentials = _credential
}

object CustomKeyMode

final case class CustomKeyMode(override val name: String,
                                    accessKey: String,
                                    secretKey: String,
                                    region: Option[String]
                                    ) extends AwsAuthMode {
  private lazy val _credential: AwsCredentials = {
    // Validate credentials synchronously here, without retry.
    // It's very unlikely to fail as it should not happen more than a few times
    // (one for the engine and for each backend using it) per Cromwell instance.
    validateCredential(AwsBasicCredentials.create(accessKey, secretKey), region)
  }

  override def credential(options: OptionLookup): AwsCredentials = _credential
}

final case class DefaultMode(override val name: String, region: Option[String]) extends AwsAuthMode {
  private lazy val _credential: AwsCredentials = {
    //
    // The ProfileCredentialsProvider will return your [default]
    // credential profile by reading from the credentials file located at
    // (~/.aws/credentials).
    //

    // Validate credentials synchronously here, without retry.
    // It's very unlikely to fail as it should not happen more than a few times
    // (one for the engine and for each backend using it) per Cromwell instance.
    validateCredential(DefaultCredentialsProvider.create.resolveCredentials(), region)
  }

  override def credential(options: OptionLookup): AwsCredentials = _credential
}


final case class AssumeRoleMode(override val name: String,
                          baseAuthName: String,
                          roleArn: String,
                          externalId: String,
                          region: Option[String]
                          ) extends AwsAuthMode {

  private lazy val _credential: AwsCredentials = {
    val requestBuilder = AssumeRoleRequest
                           .builder
                           .roleSessionName("cromwell")
                           .roleArn(roleArn)
                           .durationSeconds(3600)

    // The builder is simply mutating itself (TODO: find good ref, as v2
    // uses generated code)
    // So we can get away with a val and discard the return value
    if (! externalId.isEmpty) requestBuilder.externalId(externalId)
    val request = requestBuilder.build

    val builder = StsClient.builder
    region.foreach(str => builder.region(Region.of(str)))
    // See comment above regarding builder
    baseAuthObj match{
      case Some(auth) => builder.credentialsProvider(StaticCredentialsProvider.create(auth.credential(_ => "")))
      case _ => throw new RuntimeException(s"Base auth configuration required for assume role")
    }

    val stsCredentials = builder.build.assumeRole(request).credentials

    val sessionCredentials = AwsSessionCredentials.create(
                               stsCredentials.accessKeyId,
                               stsCredentials.secretAccessKey,
                               stsCredentials.sessionToken)

    validateCredential(sessionCredentials, region)
  }

  override def credential(options: OptionLookup): AwsCredentials = _credential

  private var baseAuthObj : Option[AwsAuthMode] = None

  def assign(baseAuth: AwsAuthMode) : Unit = {
    baseAuthObj match {
      case None => baseAuthObj = Some(baseAuth)
      case _ => throw new RuntimeException(s"Base auth object has already been assigned")
    }
  }

  // We want to allow our tests access to the value
  // of the baseAuthObj
  def baseAuthentication() : AwsAuthMode = {
    baseAuthObj match {
      case Some(o) => o
      case _ => throw new RuntimeException(s"Base auth object has not been set")
    }
  }
}

class OptionLookupException(val key: String, cause: Throwable) extends RuntimeException(key, cause)
