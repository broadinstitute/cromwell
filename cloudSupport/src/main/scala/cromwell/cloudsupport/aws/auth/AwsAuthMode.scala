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
import software.amazon.awssdk.services.sts.auth.StsAssumeRoleCredentialsProvider
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

  def provider(): AwsCredentialsProvider

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

  override def provider(): AwsCredentialsProvider = AnonymousCredentialsProvider.create
}

object CustomKeyMode

/**
 * The AwsAuthMode constructed from a 'custom_key' auths scheme.
 *
 * @param name
 * @param accessKey static AWS access key
 * @param secretKey static AWS secret key
 * @param region an optional AWS region
 */
final case class CustomKeyMode(override val name: String,
                                    accessKey: String,
                                    secretKey: String,
                                    region: Option[String]
                                    ) extends AwsAuthMode {
  private lazy val _provider: AwsCredentialsProvider = {
    // Validate credentials synchronously here, without retry.
    // It's very unlikely to fail as it should not happen more than a few times
    // (one for the engine and for each backend using it) per Cromwell instance.
    StaticCredentialsProvider.create(validateCredential(AwsBasicCredentials.create(accessKey, secretKey), region))
  }

  override def provider(): AwsCredentialsProvider = _provider
}

/**
 * The AwsAuthMode constructed from a 'default' auths scheme.
 *
 * @param name
 * @param region an optional AWS region
 */
final case class DefaultMode(override val name: String, region: Option[String]) extends AwsAuthMode {
  // The DefaultCredentialsProvider will look through a chain of standard AWS providers as
  // per the normal behaviour of aws-cli etc
  private lazy val _provider: AwsCredentialsProvider = {
    DefaultCredentialsProvider.create()
  }

  override def provider(): AwsCredentialsProvider = _provider
}

/**
 * The AwsAuthMode constructed from an 'assume_role' auths scheme.
 *
 * @param name
 * @param baseAuthName the name of another AwsAuthMode that will give us starting credentials
 * @param roleArn the ARN of the role that we want to assume
 * @param externalId an optional external id
 * @param region an optional AWS region
 */
final case class AssumeRoleMode(override val name: String,
                          baseAuthName: String,
                          roleArn: String,
                          externalId: String,
                          region: Option[String]
                          ) extends AwsAuthMode {

  private lazy val _provider: StsAssumeRoleCredentialsProvider = {
    // we need to perform operations on STS using the credentials provided from the baseAuthName
    val stsBuilder = StsClient.builder
    region.foreach(str => stsBuilder.region(Region.of(str)))
    baseAuthObj match{
      case Some(auth) => stsBuilder.credentialsProvider(auth.provider())
      case _ => throw new RuntimeException(s"Base auth configuration required for assume role")
    }

    // the STS operation we are going to perform is an assume-role to the given role ARN
    val assumeRoleBuilder = AssumeRoleRequest.builder
      .roleArn(roleArn)
      .durationSeconds(3600)
      .roleSessionName("cromwell")
    if (! externalId.isEmpty) assumeRoleBuilder.externalId(externalId)

    // the returned provider will handle refreshing the assume-role creds when needed
    StsAssumeRoleCredentialsProvider.builder
      .stsClient(stsBuilder.build())
      .refreshRequest(assumeRoleBuilder.build())
      .build()
  }

  // we use the lazily create _provider as providers that support refreshing token will sometimes
  // start a background thread to perform the refresh
  override def provider(): AwsCredentialsProvider = _provider

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
