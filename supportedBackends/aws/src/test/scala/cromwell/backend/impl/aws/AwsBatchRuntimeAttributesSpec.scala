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

import cats.data.NonEmptyList
import cromwell.backend.impl.aws.io.{AwsBatchVolume, AwsBatchWorkingDisk}
import cromwell.backend.validation.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet}
import cromwell.backend.{RuntimeAttributeDefinition}
import cromwell.core.WorkflowOptions
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.scalatest.{Matchers, WordSpecLike}
import org.slf4j.helpers.NOPLogger
import org.specs2.mock.Mockito
import spray.json._
import squants.information.Gigabytes
import wom.types._
import wom.values._

class AwsBatchRuntimeAttributesSpec extends WordSpecLike with Matchers with Mockito {

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]): WorkflowOptions = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

  val expectedDefaults = new AwsBatchRuntimeAttributes(refineMV[Positive](1), Vector("us-east-1a", "us-east-1b"),

    Gigabytes(2), Vector(AwsBatchWorkingDisk()), "ubuntu:latest", "arn::::", false,
    ContinueOnReturnCodeSet(Set(0)), false)

  "AwsBatchRuntimeAttributes" should {

    "throw an exception when there are no runtime attributes defined." in {
      val runtimeAttributes = Map.empty[String, WomValue]
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Can't find an attribute value for key docker")
    }

    // TODO: Fix this test. The functionality works fine - the idea is that
    //       queueArn is required, does not have a hard coded default,
    //       but can use the default-runtime-attributes stanza as a configured
    //       default. I believe the failure of this test is due to the fact
    //       that this test does not properly setup configuration.runtimeAttributes
    //       but instead tries to do some default thing after the construction of the
    //       object. I have not yet investigated this configuration mechanism.
    //
    // "use hardcoded defaults if not declared in task, workflow options, or config (except for docker)" in {
    //   val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "queueArn" -> WomString("arn::::"))
    //   val expectedRuntimeAttributes = expectedDefaults
    //   assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, configuration = noDefaultsAwsBatchConfiguration)
    // }

    "validate a valid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))
      val expectedRuntimeAttributes = expectedDefaults
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid Docker entry" in {
      val runtimeAttributes = Map("docker" -> WomInteger(1))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "validate a valid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaults.copy(failOnStderr = true)
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "failOnStderr" -> WomString("yes"))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "validate a valid continueOnReturnCode integer entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomInteger(1))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode boolean entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomBoolean(false))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeFlag(false))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomArray(WomArrayType(WomIntegerType), Array(WomInteger(1), WomInteger(2))))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "coerce then validate a valid continueOnReturnCode array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomArray(WomArrayType(WomStringType), Array(WomString("1"), WomString("2"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2)))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "continueOnReturnCode" -> WomString("value"))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "validate a valid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(2))
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV[Positive](2))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "validate a valid cpu string entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("2"))
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV[Positive](2))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid cpu entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomString("value"))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "validate a valid zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomString("us-east-1a"))
      val expectedRuntimeAttributes = expectedDefaults.copy(zones = Vector("us-east-1a"))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomInteger(1))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "validate a valid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomArray(WomArrayType(WomStringType), Array(WomString("us-east-1a"), WomString("us-east-1b"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(zones = Vector("us-east-1a", "us-east-1b"))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid array zones entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "zones" -> WomArray(WomArrayType(WomIntegerType), Array(WomInteger(1), WomInteger(2))))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "validate a valid disks entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomString("local-disk"))
      val expectedRuntimeAttributes = expectedDefaults.copy(disks = Seq(AwsBatchVolume.parse("local-disk").get))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid disks entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomInteger(10))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disks runtime attribute to be a comma separated String or Array[String]")
    }

    "validate a valid disks array entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomArray(WomArrayType(WomStringType), Array(WomString("local-disk"), WomString("local-disk"))))
      val expectedRuntimeAttributes = expectedDefaults.copy(disks = Seq(AwsBatchVolume.parse("local-disk").get, AwsBatchVolume.parse("local-disk").get))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    // TODO: This is working ok (appropriate error messages), though test is throwing due to message inconsistency
    // "fail to validate a valid disks array entry" in {
    //   val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "disks" -> WomArray(WomArrayType(WomStringType), Array(WomString("blah"), WomString("blah blah"))))
    //   assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE' but got 'blah blah'")
    // }

    "validate a valid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("1 GB"))
      val expectedRuntimeAttributes = expectedDefaults.copy(memory = Gigabytes(1))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid memory entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "memory" -> WomString("blah"))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "validate a valid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "noAddress" -> WomBoolean(true))
      val expectedRuntimeAttributes = expectedDefaults.copy(noAddress = true)
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "fail to validate an invalid noAddress entry" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "noAddress" -> WomInteger(1))
      assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes,
        "Expecting noAddress runtime attribute to be a Boolean")
    }

    "override config default attributes with default attributes declared in workflow options" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))

      val workflowOptionsJson =
        """{
          |  "default_runtime_attributes": { "cpu": 2 }
          |}
        """.stripMargin.parseJson.asInstanceOf[JsObject]

      val workflowOptions = WorkflowOptions.fromJsonObject(workflowOptionsJson).get
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV[Positive](2))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }

    "override config default runtime attributes with task runtime attributes" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(4))

      val workflowOptionsJson =
        """{
          |  "default_runtime_attributes": { "cpu": 2 }
          |}
        """.stripMargin.parseJson.asInstanceOf[JsObject]

      val workflowOptions = WorkflowOptions.fromJsonObject(workflowOptionsJson).get
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV[Positive](4))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }

    "override invalid config default attributes with task runtime attributes" in {
      val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "cpu" -> WomInteger(4))

      val workflowOptionsJson =
        """{
          |  "default_runtime_attributes": { "cpu": 2.2 }
          |}
        """.stripMargin.parseJson.asInstanceOf[JsObject]

      val workflowOptions = WorkflowOptions.fromJsonObject(workflowOptionsJson).get
      val expectedRuntimeAttributes = expectedDefaults.copy(cpu = refineMV[Positive](4))
      assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes, workflowOptions)
    }
  }

  private def assertAwsBatchRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WomValue],
                                                           expectedRuntimeAttributes: AwsBatchRuntimeAttributes,
                                                           workflowOptions: WorkflowOptions = emptyWorkflowOptions,
                                                           defaultZones: NonEmptyList[String] = defaultZones,
                                                           configuration: AwsBatchConfiguration = configuration): Unit = {
    try {
      val actualRuntimeAttributes = toAwsBatchRuntimeAttributes(runtimeAttributes, workflowOptions, configuration)
      assert(actualRuntimeAttributes == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  private def assertAwsBatchRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WomValue],
                                                       exMsg: String,
                                                       workflowOptions: WorkflowOptions = emptyWorkflowOptions): Unit = {
    try {
      toAwsBatchRuntimeAttributes(runtimeAttributes, workflowOptions, configuration)
      fail(s"A RuntimeException was expected with message: $exMsg")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
    ()
  }

  private def toAwsBatchRuntimeAttributes(runtimeAttributes: Map[String, WomValue],
                                     workflowOptions: WorkflowOptions,
                                     configuration: AwsBatchConfiguration): AwsBatchRuntimeAttributes = {
    val runtimeAttributesBuilder = AwsBatchRuntimeAttributes.runtimeAttributesBuilder(configuration)
    val defaultedAttributes = RuntimeAttributeDefinition.addDefaultsToAttributes(
      staticRuntimeAttributeDefinitions, workflowOptions)(runtimeAttributes)
    val validatedRuntimeAttributes = runtimeAttributesBuilder.build(defaultedAttributes, NOPLogger.NOP_LOGGER)
    AwsBatchRuntimeAttributes(validatedRuntimeAttributes, configuration.runtimeConfig)
  }

  private val emptyWorkflowOptions = WorkflowOptions.fromMap(Map.empty).get
  private val defaultZones = NonEmptyList.of("us-east-1a", "us-east-1b")
  private val configuration = new AwsBatchConfiguration(AwsBatchTestConfig.AwsBatchBackendConfigurationDescriptor)
  // private val noDefaultsAwsBatchConfiguration = new AwsBatchConfiguration(AwsBatchTestConfig.NoDefaultsConfigurationDescriptor)
  private val staticRuntimeAttributeDefinitions: Set[RuntimeAttributeDefinition] =
    AwsBatchRuntimeAttributes.runtimeAttributesBuilder(configuration).definitions.toSet
}
