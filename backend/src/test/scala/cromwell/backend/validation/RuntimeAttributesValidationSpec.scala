package cromwell.backend.validation

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.TestConfig
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import wom.RuntimeAttributesKeys
import wom.types._
import wom.values._

class RuntimeAttributesValidationSpec extends WordSpecLike with Matchers with BeforeAndAfterAll {

  val mockBackendRuntimeConfig = TestConfig.allRuntimeAttrsConfig

  "RuntimeAttributesValidation" should {
    "return success when tries to validate a valid Docker entry" in {
      val dockerValue = Some(WomString("someImage"))
      val result = RuntimeAttributesValidation.validateDocker(dockerValue,
        "Failed to get Docker mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x.get == "someImage")
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success (based on defined HoF) when tries to validate a docker entry but it does not contain a value" in {
      val dockerValue = None
      val result = RuntimeAttributesValidation.validateDocker(dockerValue, None.validNel)
      result match {
        case Valid(x) => assert(x.isEmpty)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure (based on defined HoF) when tries to validate a docker entry but it does not contain a value" in {
      val dockerValue = None
      val result = RuntimeAttributesValidation.validateDocker(dockerValue,
        "Failed to get Docker mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Failed to get Docker mandatory key from runtime attributes")
      }
    }

    "return failure when there is an invalid docker runtime attribute defined" in {
      val dockerValue = Some(WomInteger(1))
      val result = RuntimeAttributesValidation.validateDocker(dockerValue,
        "Failed to get Docker mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting docker runtime attribute to be a String")
      }
    }

    "return success when tries to validate a failOnStderr boolean entry" in {
      val failOnStderrValue = Some(WomBoolean(true))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        "Failed to get failOnStderr mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a failOnStderr 'true' string entry" in {
      val failOnStderrValue = Some(WomString("true"))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        "Failed to get failOnStderr mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a failOnStderr 'false' string entry" in {
      val failOnStderrValue = Some(WomString("false"))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        "Failed to get failOnStderr mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(!x)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid failOnStderr runtime attribute defined" in {
      val failOnStderrValue = Some(WomInteger(1))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        "Failed to get failOnStderr mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
      }
    }

    "return success (based on defined HoF) when tries to validate a failOnStderr entry but it does not contain a value" in {
      val failOnStderrValue = None
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue, true.validNel)
      result match {
        case Valid(x) => assert(x)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode boolean entry" in {
      val continueOnReturnCodeValue = Some(WomBoolean(true))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeFlag(true))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode 'true' string entry" in {
      val continueOnReturnCodeValue = Some(WomString("true"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeFlag(true))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode 'false' string entry" in {
      val continueOnReturnCodeValue = Some(WomString("false"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeFlag(false))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode int entry" in {
      val continueOnReturnCodeValue = Some(WomInteger(12))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeSet(Set(12)))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid continueOnReturnCode runtime attribute defined" in {
      val continueOnReturnCodeValue = Some(WomString("yes"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) =>
          assert(e.head == "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
      }
    }

    "return success when there is a valid integer array in continueOnReturnCode runtime attribute" in {
      val continueOnReturnCodeValue = Some(WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(1), WomInteger(2))))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeSet(Set(1, 2)))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid array in continueOnReturnCode runtime attribute" in {
      val continueOnReturnCodeValue = Some(WomArray(WomArrayType(WomStringType), Seq(WomString("one"), WomString("two"))))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
      }
    }

    "return success (based on defined HoF) when tries to validate a continueOnReturnCode entry but it does not contain a value" in {
      val continueOnReturnCodeValue = None
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue, ContinueOnReturnCodeFlag(false).validNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeFlag(false))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a valid Integer memory entry" in {
      val expectedGb = 1
      val memoryValue = Some(WomInteger(1 << 30))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x.amount == expectedGb)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when tries to validate an invalid Integer memory entry" in {
      val memoryValue = Some(WomInteger(-1))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting memory runtime attribute value greater than 0 but got -1")
      }
    }

    "return success when tries to validate a valid String memory entry" in {
      val expectedGb = 2
      val memoryValue = Some(WomString("2 GB"))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x.amount == expectedGb)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when tries to validate an invalid size in String memory entry" in {
      val memoryValue = Some(WomString("0 GB"))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting memory runtime attribute value greater than 0 but got 0.0")
      }
    }

    "return failure when tries to validate an invalid String memory entry" in {
      val memoryValue = Some(WomString("value"))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting memory runtime attribute to be an Integer or String with format '8 GB'. Exception: value should be of the form 'X Unit' where X is a number, e.g. 8 GB")
      }
    }

    "return failure when tries to validate an invalid memory entry" in {
      val memoryValue = Some(WomBoolean(true))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting memory runtime attribute to be an Integer or String with format '8 GB'. Exception: Not supported WDL type value")
      }
    }

    "return failure when tries to validate a non-provided memory entry" in {
      val memoryValue = None
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Failed to get memory mandatory key from runtime attributes")
      }
    }

    "return success when tries to validate a valid cpu entry" in {
      val cpuValue = Some(WomInteger(1))
      val result = RuntimeAttributesValidation.validateCpu(cpuValue,
        "Failed to get cpu mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x.value == 1)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when tries to validate an invalid cpu entry" in {
      val cpuValue = Some(WomInteger(-1))
      val result = RuntimeAttributesValidation.validateCpu(cpuValue,
        "Failed to get cpu mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting cpu runtime attribute value greater than 0")
      }
    }

    "return failure when tries to validate a non-provided cpu entry" in {
      val cpuValue = None
      val result = RuntimeAttributesValidation.validateCpu(cpuValue,
        "Failed to get cpu mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Failed to get cpu mandatory key from runtime attributes")
      }
    }

    "return default values as WdlValues when they can be coerced into expected WdlTypes" in {
      val optionalConfig = Option(TestConfig.allRuntimeAttrsConfig)

      val defaultVals = Map(
        "cpu" -> CpuValidation.configDefaultWomValue(optionalConfig).get,
        "failOnStderr" -> FailOnStderrValidation.configDefaultWdlValue(optionalConfig).get,
        "continueOnReturnCode" -> ContinueOnReturnCodeValidation.configDefaultWdlValue(optionalConfig).get
      )

      val expectedDefaultVals = Map(
        "cpu" -> WomInteger(1),
        "failOnStderr" -> WomBoolean(false),
        "continueOnReturnCode" -> WomInteger(0)
      )

      defaultVals shouldBe expectedDefaultVals
    }

    "return default values as BadDefaultAttribute when they can't be coerced to expected WdlTypes" in {
     val optionalInvalidAttrsConfig = Option(ConfigFactory.parseString(
       """
         |cpu = 1.4
         |failOnStderr = "notReal"
         |continueOnReturnCode = 0
       """.stripMargin))

     val defaultVals = Map(
       "cpu" -> CpuValidation.configDefaultWomValue(optionalInvalidAttrsConfig).get,
       "failOnStderr" -> FailOnStderrValidation.configDefaultWdlValue(optionalInvalidAttrsConfig).get,
       "continueOnReturnCode" -> ContinueOnReturnCodeValidation.configDefaultWdlValue(optionalInvalidAttrsConfig).get
     )

     val expectedDefaultVals = Map(
       "cpu" -> BadDefaultAttribute(WomString("1.4")),
       "failOnStderr" -> BadDefaultAttribute(WomString("notReal")),
       "continueOnReturnCode" -> WomInteger(0)
     )

     defaultVals shouldBe expectedDefaultVals
    }

    "should parse memory successfully" in {
      val backendConfigTemplate: String =
        s"""
           |  default-runtime-attributes {
           |     memory: "2 GB"
           |     memoryMin: "0.3 GB"
           |     memoryMax: "0.4 GB"
           |  }
           |""".stripMargin

      val backendConfig: Config = ConfigFactory.parseString(backendConfigTemplate).getConfig("default-runtime-attributes")

      val memoryVal = MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryKey, Some(backendConfig))
      val memoryMinVal = MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryMinKey, Some(backendConfig))
      val memoryMaxVal = MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryMaxKey, Some(backendConfig))
      MemoryValidation.withDefaultMemory(RuntimeAttributesKeys.MemoryKey, memoryVal.get).runtimeAttributeDefinition.factoryDefault shouldBe Some((WomLong(2147483648L)))
      MemoryValidation.withDefaultMemory(RuntimeAttributesKeys.MemoryMinKey, memoryMinVal.get).runtimeAttributeDefinition.factoryDefault shouldBe Some((WomLong(322122547L)))
      MemoryValidation.withDefaultMemory(RuntimeAttributesKeys.MemoryMaxKey, memoryMaxVal.get).runtimeAttributeDefinition.factoryDefault shouldBe Some((WomLong(429496729L)))
    }

    "shouldn't throw up if the value for a default-runtime-attribute key cannot be coerced into an expected WomType" in {
      val backendConfigTemplate: String =
        s"""
           |  default-runtime-attributes {
           |     memory: "blahblah"
           |  }
           |""".stripMargin

      val backendConfig: Config = ConfigFactory.parseString(backendConfigTemplate).getConfig("default-runtime-attributes")

      val memoryVal = MemoryValidation.configDefaultString(RuntimeAttributesKeys.MemoryKey, Some(backendConfig))
      MemoryValidation.withDefaultMemory(RuntimeAttributesKeys.MemoryKey, memoryVal.get).runtimeAttributeDefinition.factoryDefault shouldBe Some(BadDefaultAttribute(WomString("blahblah")))
    }

    "should be able to coerce a list of return codes into an WdlArray" in {
      val optinalBackendConfig = Option(ConfigFactory.parseString(
        s"""
           |continueOnReturnCode = [0,1,2]
           |""".stripMargin))

      ContinueOnReturnCodeValidation.configDefaultWdlValue(optinalBackendConfig).get shouldBe WomArray(WomArrayType(WomIntegerType), Array(WomInteger(0), WomInteger(1), WomInteger(2)))
    }

    "return failure when tries to validate an invalid maxRetries entry" in {
      val maxRetries = Option(WomInteger(-1))
      val result = RuntimeAttributesValidation.validateMaxRetries(maxRetries,
        "Failed to get maxRetries key from runtime attributes".invalidNel)
      result match {
        case Valid(_) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting maxRetries runtime attribute value greater than or equal to 0")
      }
    }
  }
}
