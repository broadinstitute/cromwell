package cromwell.backend.validation

import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import wdl4s.types.{WdlIntegerType, WdlStringType, WdlArrayType}
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString}

import scalaz.Scalaz._

class RuntimeAttributesValidationSpec extends WordSpecLike with Matchers with BeforeAndAfterAll {

  "RuntimeAttributesValidation" should {
    "return success when tries to validate a valid Docker entry" in {
      val dockerValue = Some(WdlString("someImage"))
      val result = RuntimeAttributesValidation.validateDocker(dockerValue,
        () => "Failed to get Docker mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(x.get == "someImage")
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success (based on defined HoF) when tries to validate a docker entry but it does not contain a value" in {
      val dockerValue = None
      val result = RuntimeAttributesValidation.validateDocker(dockerValue, () => None.successNel)
      result match {
        case scalaz.Success(x) => assert(!x.isDefined)
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure (based on defined HoF) when tries to validate a docker entry but it does not contain a value" in {
      val dockerValue = None
      val result = RuntimeAttributesValidation.validateDocker(dockerValue,
        () => "Failed to get Docker mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => fail("A failure was expected.")
        case scalaz.Failure(e) => assert(e.head == "Failed to get Docker mandatory key from runtime attributes")
      }
    }

    "return failure when there is an invalid docker runtime attribute defined" in {
      val dockerValue = Some(WdlInteger(1))
      val result = RuntimeAttributesValidation.validateDocker(dockerValue,
        () => "Failed to get Docker mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => fail("A failure was expected.")
        case scalaz.Failure(e) => assert(e.head == "Expecting docker runtime attribute to be a String")
      }
    }

    "return success when tries to validate a failOnStderr boolean entry" in {
      val failOnStderrValue = Some(WdlBoolean(true))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        () => "Failed to get failOnStderr mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(x)
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a failOnStderr 'true' string entry" in {
      val failOnStderrValue = Some(WdlString("true"))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        () => "Failed to get failOnStderr mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(x)
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a failOnStderr 'false' string entry" in {
      val failOnStderrValue = Some(WdlString("false"))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        () => "Failed to get failOnStderr mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(!x)
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid failOnStderr runtime attribute defined" in {
      val failOnStderrValue = Some(WdlInteger(1))
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue,
        () => "Failed to get failOnStderr mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => fail("A failure was expected.")
        case scalaz.Failure(e) => assert(e.head == "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
      }
    }

    "return success (based on defined HoF) when tries to validate a failOnStderr entry but it does not contain a value" in {
      val failOnStderrValue = None
      val result = RuntimeAttributesValidation.validateFailOnStderr(failOnStderrValue, () => true.successNel)
      result match {
        case scalaz.Success(x) => assert(x)
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode boolean entry" in {
      val continueOnReturnCodeValue = Some(WdlBoolean(true))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        () => "Failed to get continueOnReturnCode mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(x == ContinueOnReturnCodeFlag(true))
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode 'true' string entry" in {
      val continueOnReturnCodeValue = Some(WdlString("true"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        () => "Failed to get continueOnReturnCode mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(x == ContinueOnReturnCodeFlag(true))
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode 'false' string entry" in {
      val continueOnReturnCodeValue = Some(WdlString("false"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        () => "Failed to get continueOnReturnCode mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(x == ContinueOnReturnCodeFlag(false))
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return success when tries to validate a continueOnReturnCode int entry" in {
      val continueOnReturnCodeValue = Some(WdlInteger(12))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        () => "Failed to get continueOnReturnCode mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(x == ContinueOnReturnCodeSet(Set(12)))
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid continueOnReturnCode runtime attribute defined" in {
      val continueOnReturnCodeValue = Some(WdlString("yes"))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        () => "Failed to get continueOnReturnCode mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => fail("A failure was expected.")
        case scalaz.Failure(e) =>
          assert(e.head == "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
      }
    }

    "return success when there is a valid integer array in continueOnReturnCode runtime attribute" in {
      val continueOnReturnCodeValue = Some(WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(1), WdlInteger(2))))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        () => "Failed to get continueOnReturnCode mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => assert(x == ContinueOnReturnCodeSet(Set(1, 2)))
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }

    "return failure when there is an invalid array in continueOnReturnCode runtime attribute" in {
      val continueOnReturnCodeValue = Some(WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("one"), WdlString("two"))))
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue,
        () => "Failed to get continueOnReturnCode mandatory key from runtime attributes".failureNel)
      result match {
        case scalaz.Success(x) => fail("A failure was expected.")
        case scalaz.Failure(e) => assert(e.head == "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
      }
    }

    "return success (based on defined HoF) when tries to validate a continueOnReturnCode entry but it does not contain a value" in {
      val continueOnReturnCodeValue = None
      val result = RuntimeAttributesValidation.validateContinueOnReturnCode(continueOnReturnCodeValue, () => ContinueOnReturnCodeFlag(false).successNel)
      result match {
        case scalaz.Success(x) => assert(x == ContinueOnReturnCodeFlag(false))
        case scalaz.Failure(e) => fail(e.toList.mkString(" "))
      }
    }
  }
}
