package cromwell.backend.validation

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString}

class RuntimeAttributesValidationSpec extends WordSpecLike with Matchers with BeforeAndAfterAll {

  "RuntimeAttributesValidation" should {
    "return success when tries to validate a valid Docker entry" in {
      val dockerValue = Some(WdlString("someImage"))
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
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Failed to get Docker mandatory key from runtime attributes")
      }
    }

    "return failure when there is an invalid docker runtime attribute defined" in {
      val dockerValue = Some(WdlInteger(1))
      val result = RuntimeAttributesValidation.validateDocker(dockerValue,
        "Failed to get Docker mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting docker runtime attribute to be a String")
      }
    }

    "return success when tries to validate a failOnStderr boolean entry" in {
      val failOnStderrValue = Some(WdlBoolean(true))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        "Failed to get failOnStderr mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a failOnStderr 'true' string entry" in {
      val failOnStderrValue = Some(WdlString("true"))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        "Failed to get failOnStderr mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a failOnStderr 'false' string entry" in {
      val failOnStderrValue = Some(WdlString("false"))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        "Failed to get failOnStderr mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(!x)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid failOnStderr runtime attribute defined" in {
      val failOnStderrValue = Some(WdlInteger(1))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        "Failed to get failOnStderr mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
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
      val continueOnReturnCodeValue = Some(WdlBoolean(true))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeFlag(true))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode 'true' string entry" in {
      val continueOnReturnCodeValue = Some(WdlString("true"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeFlag(true))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode 'false' string entry" in {
      val continueOnReturnCodeValue = Some(WdlString("false"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeFlag(false))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode int entry" in {
      val continueOnReturnCodeValue = Some(WdlInteger(12))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeSet(Set(12)))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid continueOnReturnCode runtime attribute defined" in {
      val continueOnReturnCodeValue = Some(WdlString("yes"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) =>
          assert(e.head == "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
      }
    }

    "return success when there is a valid integer array in continueOnReturnCode runtime attribute" in {
      val continueOnReturnCodeValue = Some(WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(1), WdlInteger(2))))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == ContinueOnReturnCodeSet(Set(1, 2)))
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid array in continueOnReturnCode runtime attribute" in {
      val continueOnReturnCodeValue = Some(WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("one"), WdlString("two"))))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        "Failed to get continueOnReturnCode mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
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
      val memoryValue = Some(WdlInteger(1000000000))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x.amount == expectedGb)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when tries to validate an invalid Integer memory entry" in {
      val memoryValue = Some(WdlInteger(-1))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting memory runtime attribute value greater than 0 but got -1")
      }
    }

    "return success when tries to validate a valid String memory entry" in {
      val expectedGb = 2
      val memoryValue = Some(WdlString("2 GB"))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x.amount == expectedGb)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when tries to validate an invalid size in String memory entry" in {
      val memoryValue = Some(WdlString("0 GB"))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting memory runtime attribute value greater than 0 but got 0.0")
      }
    }

    "return failure when tries to validate an invalid String memory entry" in {
      val memoryValue = Some(WdlString("value"))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting memory runtime attribute to be an Integer or String with format '8 GB'. Exception: value should be of the form 'X Unit' where X is a number, e.g. 8 GB")
      }
    }

    "return failure when tries to validate an invalid memory entry" in {
      val memoryValue = Some(WdlBoolean(true))
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting memory runtime attribute to be an Integer or String with format '8 GB'. Exception: Not supported WDL type value")
      }
    }

    "return failure when tries to validate a non-provided memory entry" in {
      val memoryValue = None
      val result = RuntimeAttributesValidation.validateMemory(memoryValue,
        "Failed to get memory mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Failed to get memory mandatory key from runtime attributes")
      }
    }

    "return success when tries to validate a valid cpu entry" in {
      val cpuValue = Some(WdlInteger(1))
      val result = RuntimeAttributesValidation.validateCpu(cpuValue,
        "Failed to get cpu mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => assert(x == 1)
        case Invalid(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when tries to validate an invalid cpu entry" in {
      val cpuValue = Some(WdlInteger(-1))
      val result = RuntimeAttributesValidation.validateCpu(cpuValue,
        "Failed to get cpu mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Expecting cpu runtime attribute value greater than 0")
      }
    }

    "return failure when tries to validate a non-provided cpu entry" in {
      val cpuValue = None
      val result = RuntimeAttributesValidation.validateMemory(cpuValue,
        "Failed to get cpu mandatory key from runtime attributes".invalidNel)
      result match {
        case Valid(x) => fail("A failure was expected.")
        case Invalid(e) => assert(e.head == "Failed to get cpu mandatory key from runtime attributes")
      }
    }
  }
}
