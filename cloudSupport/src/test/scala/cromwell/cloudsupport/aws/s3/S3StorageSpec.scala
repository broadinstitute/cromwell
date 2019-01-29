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
package cromwell.cloudsupport.aws.s3
import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import org.scalatest.{FlatSpec, Matchers, Tag}
import software.amazon.awssdk.regions.Region

class S3StorageSpec extends FlatSpec with Matchers {

  behavior of "S3Storage"

  it should "build the default cloud storage configuration" taggedAs S3StorageSpecUtils.AwsTest in {
    val configuration = S3Storage.DefaultConfiguration
    configuration.accelerateModeEnabled should be(false)
    configuration.dualstackEnabled should be(false)
    configuration.pathStyleAccessEnabled should be(false)
  }

  it should "build s3 storage" taggedAs S3StorageSpecUtils.AwsTest in {
    val configuration = S3Storage.s3Configuration(dualstackEnabled = true)
    configuration.accelerateModeEnabled should be(false)
    configuration.dualstackEnabled should be(true)
    configuration.pathStyleAccessEnabled should be(false)
  }

  it should "build s3 client with credentials" taggedAs S3StorageSpecUtils.AwsTest in {
    S3Storage.s3Client(
      AnonymousCredentialsProvider.create.resolveCredentials(),
      Option(Region.US_EAST_1)
    )
  }

}

object S3StorageSpecUtils {
  val AwsTest = Tag("AwsTest")
}
