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

import cats.implicits._
import com.typesafe.config.{ConfigException, ConfigFactory}
import cromwell.cloudsupport.aws.AwsConfiguration.AwsConfigurationException
import cromwell.cloudsupport.aws.auth.{AssumeRoleMode,CustomKeyMode,DefaultMode}
import org.scalatest.{FlatSpec, Matchers}


class AwsConfigurationSpec extends FlatSpec with Matchers {

  behavior of "AwsConfiguration"

  it should "parse all manner of well-formed auths" in {
    val righteousAwsConfig =
      s"""
        |aws {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "default"
        |      scheme = "default"
        |    },
        |    {
        |      name = "custom-keys"
        |      scheme = "custom_keys"
        |      access-key = "access_key_id"
        |      secret-key = "secret_key"
        |    },
        |    {
        |      name = "assume-role-based-on-another-with-external"
        |      scheme = "assume_role"
        |      base-auth = "default"
        |      role-arn = "my-role-arn"
        |      external-id = "my-external-id"
        |    },
        |    {
        |      name = "assume-role-based-on-another"
        |      scheme = "assume_role"
        |      base-auth = "default"
        |      role-arn = "my-role-arn"
        |    }
        |  ]
        |
        |  region = "region"
        |}
        |
      """.stripMargin

    val conf = AwsConfiguration(ConfigFactory.parseString(righteousAwsConfig))

    conf.applicationName shouldBe "cromwell"
    conf.region shouldBe "region"
    conf.authsByName should have size 4

    val auths = conf.authsByName.values

    val default = (auths collectFirst { case a: DefaultMode => a }).get
    default.name shouldBe "default"

    val customKey = (auths collectFirst { case a: CustomKeyMode => a }).get
    customKey.name shouldBe "custom-keys"
    customKey.accessKey shouldBe "access_key_id"
    customKey.secretKey shouldBe "secret_key"

    val assumeRoleWithId = (auths collectFirst { case a: AssumeRoleMode => a }).get
    assumeRoleWithId.name shouldBe "assume-role-based-on-another-with-external"
    assumeRoleWithId.baseAuthName shouldBe "default"
    assumeRoleWithId.baseAuthentication.name shouldBe "default"
    assumeRoleWithId.roleArn shouldBe "my-role-arn"
    assumeRoleWithId.externalId shouldBe "my-external-id"

    val assumeRole = (auths.takeRight(1) collectFirst { case a: AssumeRoleMode => a }).get
    assumeRole.name shouldBe "assume-role-based-on-another"
    assumeRole.baseAuthName shouldBe "default"
    assumeRole.baseAuthentication.name shouldBe "default"
    assumeRole.roleArn shouldBe "my-role-arn"
    assumeRole.externalId shouldBe ""
  }

  it should "default region to us-east-1" in {
    val config =
      """|aws {
         |  application-name = "cromwell"
         |
         |  auths = [
         |    {
         |      name = "name-default"
         |      scheme = "default"
         |    }
         |  ]
         |}
         |""".stripMargin

    val conf = AwsConfiguration(ConfigFactory.parseString(config))
    conf.region shouldBe "us-east-1"
  }

  it should "return a known auth" in {
    val config =
      """|aws {
         |  application-name = "cromwell"
         |
         |  auths = [
         |    {
         |      name = "name-default"
         |      scheme = "default"
         |    }
         |  ]
         |}
         |""".stripMargin

    val conf = AwsConfiguration(ConfigFactory.parseString(config))
    conf.auth("name-default").map(_.name) should be("name-default".valid)
  }

  it should "not return an unknown auth" in {
    val config =
      """|aws {
         |  application-name = "cromwell"
         |
         |  auths = [
         |    {
         |      name = "name-default"
         |      scheme = "default"
         |    }
         |  ]
         |}
         |""".stripMargin

    val conf = AwsConfiguration(ConfigFactory.parseString(config))
    conf.auth("name-botched") should be(
      "`aws` configuration stanza does not contain an auth named 'name-botched'.  Known auth names: name-default"
        .invalidNel)
  }

  it should "not parse a configuration stanza without applicationName" in {
    val applessAwsConfig =
      """
        |aws {
        |  auths = [
        |    {
        |      name = "name-default"
        |      scheme = "default"
        |    }
        |  ]
        |}
      """.stripMargin

    the[AwsConfigurationException] thrownBy {
      AwsConfiguration(ConfigFactory.parseString(applessAwsConfig))
    } should have message "AWS configuration:\nNo configuration setting found for key 'application-name'"
  }

  it should "not parse a configuration stanza without service account credentials" in {
    val noServiceAccountCredentials =
      """
        |aws {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "service-account"
        |      scheme = "custom_keys"
        |    }
        |  ]
        |}
      """.stripMargin

    the[AwsConfigurationException] thrownBy {
      AwsConfiguration(ConfigFactory.parseString(noServiceAccountCredentials))
    } should have message "AWS configuration:\n" +
      "Access key and/or secret key missing for service account \"service-account\". See reference.conf under the " +
      "aws.auth, custom key section for details of required configuration."
  }

  it should "not parse a configuration stanza with an unsupported authentication scheme" in {
    val unsupported =
      """
        |aws {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "unsupported-auth"
        |      scheme = "not supported"
        |    }
        |  ]
        |}
      """.stripMargin

    the[AwsConfigurationException] thrownBy {
      AwsConfiguration(ConfigFactory.parseString(unsupported))
    } should have message "AWS configuration:\nUnsupported authentication scheme: not supported"
  }

  it should "not parse a configuration stanza without a schema" in {
    val schemeless =
      """
        |aws {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "scheme-unspecified"
        |    }
        |  ]
        |}
      """.stripMargin

    the[ConfigException.Missing] thrownBy {
      AwsConfiguration(ConfigFactory.parseString(schemeless))
    } should have message "No configuration setting found for key 'scheme'"
  }

  it should "not parse a configuration stanza without an auth name" in {
    val nameless =
      """
        |aws {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      scheme = "default"
        |    }
        |  ]
        |}
      """.stripMargin

    the[ConfigException.Missing] thrownBy {
      AwsConfiguration(ConfigFactory.parseString(nameless))
    } should have message "No configuration setting found for key 'name'"
  }

  it should "not parse a configuration stanza with a bad access-key in custom keys mode" in {
    // The various AwsAuthModes actually don't complain about spurious keys in their
    // configurations as long as all the keys they do care about are present.  That's not
    // necessarily ideal behavior.
    val badKeyInRefreshTokenMode =
      """
        |aws {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "name-refresh"
        |      scheme = "custom_keys"
        |      access-key-botched-key = "secret_id"
        |      secret-key = "secret_secret"
        |    }
        |  ]
        |}
      """.stripMargin

    the[AwsConfigurationException] thrownBy {
      AwsConfiguration(ConfigFactory.parseString(badKeyInRefreshTokenMode))
    } should have message "AWS configuration:\nAccess key and/or secret key missing for service account \"name-refresh\". See reference.conf under the aws.auth, custom key section for details of required configuration."
  }

  it should "not parse a configuration stanza without a role-arn in assume-role mode" in {
    val badKeyInUserMode =
      """
        |aws {
        |  application-name = "cromwell"
        |
        |  auths = [
        |    {
        |      name = "name-user"
        |      scheme = "assume_role"
        |      role-arn-botched = "my-role-arn"
        |      base-auth = "default"
        |    }
        |  ]
        |}
      """.stripMargin

    the[AwsConfigurationException] thrownBy {
      AwsConfiguration(ConfigFactory.parseString(badKeyInUserMode))
    } should have message "AWS configuration:\nNo configuration setting found for key 'role-arn'"
  }

  it should "not parse a configuration stanza with a duplicate auth name" in {
    val duplicateAuthName =
      """|aws {
         |  application-name = "cromwell"
         |
         |  auths = [
         |    {
         |      name = "name-default"
         |      scheme = "default"
         |    }
         |    {
         |      name = "name-default"
         |      scheme = "default"
         |    }
         |  ]
         |}
         |""".stripMargin

    the[AwsConfigurationException] thrownBy {
      AwsConfiguration(ConfigFactory.parseString(duplicateAuthName))
    } should have message "AWS configuration:\nDuplicate auth names: name-default"
  }
}
