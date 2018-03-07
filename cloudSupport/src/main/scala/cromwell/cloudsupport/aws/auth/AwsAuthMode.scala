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

import java.io.{ByteArrayInputStream, FileNotFoundException}

import com.amazonaws.{AmazonClientException, AmazonServiceException}
import com.amazonaws.auth.{BasicAWSCredentials, AnonymousAWSCredentials,
                           AWSCredentials, AWSStaticCredentialsProvider,
                           DefaultAWSCredentialsProviderChain, BasicSessionCredentials}
import com.amazonaws.regions.{Region, Regions}
import com.amazonaws.services.securitytoken.{AWSSecurityTokenService, AWSSecurityTokenServiceClientBuilder}
import com.amazonaws.services.securitytoken.model.{GetCallerIdentityRequest, GetCallerIdentityResult, AssumeRoleRequest}

import cromwell.cloudsupport.aws.auth.AwsAuthMode.OptionLookup

import org.slf4j.LoggerFactory
import com.google.api.client.json.jackson2.JacksonFactory
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

  def credential(options: OptionLookup): AWSCredentials

  /**
    * Enables swapping out credential validation for various testing purposes ONLY.
    *
    * All traits in this file are sealed, all classes final, meaning things
    * like Mockito or other java/scala overrides cannot work.
    */
   private[auth] var credentialValidation: ((AWSCredentials, String) => Unit) =
     (credentials: AWSCredentials, region: String) => {
       AWSSecurityTokenServiceClientBuilder
          .standard
          .withCredentials(new AWSStaticCredentialsProvider(credentials))
          .withRegion(region)
          .build
          .getCallerIdentity(new GetCallerIdentityRequest())
       ()
     }

  protected def validateCredential(credential: AWSCredentials, region: String) = {
    Try(credentialValidation(credential, region)) match {
      case Failure(ex) => throw new RuntimeException(s"Credentials are invalid: ${ex.getMessage}", ex)
      case Success(_) => credential
    }
  }
}

case object MockAuthMode extends AwsAuthMode {
  override val name = "no_auth"

  val _credential = new AnonymousAWSCredentials

  override def credential(options: OptionLookup): AWSCredentials = _credential
}

object CustomKeyMode

final case class CustomKeyMode(override val name: String,
                                    accessKey: String,
                                    secretKey: String,
                                    region: String
                                    ) extends AwsAuthMode {
  private lazy val _credential: AWSCredentials = {
    // Validate credentials synchronously here, without retry.
    // It's very unlikely to fail as it should not happen more than a few times
    // (one for the engine and for each backend using it) per Cromwell instance.
    validateCredential(new BasicAWSCredentials(accessKey, secretKey), region)
  }

  override def credential(options: OptionLookup): AWSCredentials = _credential
}

final case class DefaultMode(override val name: String, region: String) extends AwsAuthMode {
  private lazy val _credential: AWSCredentials = {
    //
    // The ProfileCredentialsProvider will return your [default]
    // credential profile by reading from the credentials file located at
    // (~/.aws/credentials).
    //

    // Validate credentials synchronously here, without retry.
    // It's very unlikely to fail as it should not happen more than a few times
    // (one for the engine and for each backend using it) per Cromwell instance.
    validateCredential(DefaultAWSCredentialsProviderChain.getInstance.getCredentials, region)
  }

  override def credential(options: OptionLookup): AWSCredentials = _credential
}


final case class AssumeRoleMode(override val name: String,
                          baseAuth: Option[AwsAuthMode],
                          roleArn: String,
                          externalId: String,
                          region: String
                          ) extends AwsAuthMode {

  private lazy val _credential: AWSCredentials = {
    val request = new AssumeRoleRequest()
                    .withRoleSessionName("cromwell")
                    .withRoleArn(roleArn)
                    .withDurationSeconds(3600)

    // The builder is simply mutating itself (ref:
    // https://github.com/aws/aws-sdk-java/blob/4734de6fb0f80fe5768a6587aad3b9d0eaec388f/aws-java-sdk-core/src/main/java/com/amazonaws/client/builder/AwsClientBuilder.java#L395
    // So we can get away with a val and discard the return value
    if (externalId.isEmpty) request.withExternalId(externalId)

    val builder = AWSSecurityTokenServiceClientBuilder.standard.withRegion(region)
    // See comment above regarding builder
    baseAuth match{
      case Some(auth) => builder.withCredentials(new AWSStaticCredentialsProvider(auth.credential(_ => "")))
      case None => ()
    }

    val stsCredentials = builder.build.assumeRole(request).getCredentials

    val sessionCredentials = new BasicSessionCredentials(
                               stsCredentials.getAccessKeyId,
                               stsCredentials.getSecretAccessKey,
                               stsCredentials.getSessionToken)

    validateCredential(new AWSStaticCredentialsProvider(sessionCredentials).getCredentials, region)
  }

  override def credential(options: OptionLookup): AWSCredentials = _credential
}

class OptionLookupException(val key: String, cause: Throwable) extends RuntimeException(key, cause)
