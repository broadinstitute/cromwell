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

import java.io.StringWriter
import java.security.KeyPairGenerator
import java.util.Base64

import com.google.api.client.http.{HttpHeaders, HttpResponseException}
import org.scalatest.Assertions.{cancel}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

import scala.util.{Failure, Try}
import scala.language.postfixOps

class AwsAuthModeSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "AwsAuthMode"

  // TODO: Nearly everything in AwsAuthMode is going to end up validating
  //       credentials against the live service. Not much left to test
}

object AwsAuthModeSpec {
  // TODO: determine if and when this might be called
  lazy val credentialsContents: String = {
    toJson(
      "type" -> "keys",
      "access_key" -> "access_key_id",
      "secret_key" -> "secret_key"
    )
  }

  private def toJson(contents: (String, String)*): String = {
    // Generator doesn't matter as long as it generates JSON.
    // Using `jsonFactory` to get an extra line hit of coverage.
    val factory = AwsAuthMode.jsonFactory
    val writer = new StringWriter()
    val generator = factory.createJsonGenerator(writer)
    generator.enablePrettyPrint()
    generator.writeStartObject()
    contents foreach {
      case (key, value) =>
        generator.writeFieldName(key)
        generator.writeString(value)
    }
    generator.writeEndObject()
    generator.close()
    writer.toString
  }

}
