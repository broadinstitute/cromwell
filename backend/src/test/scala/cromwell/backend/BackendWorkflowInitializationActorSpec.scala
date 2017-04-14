package cromwell.backend

import akka.actor.ActorRef
import akka.testkit.TestActorRef
import com.typesafe.config.ConfigFactory
import cromwell.backend.validation.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet, ContinueOnReturnCodeValidation, RuntimeAttributesKeys}
import cromwell.core.{TestKitSuite, WorkflowOptions}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}
import wdl4s.types._
import wdl4s.values.{WdlArray, WdlBoolean, WdlFloat, WdlInteger, WdlString, WdlValue}
import wdl4s.{TaskCall, WdlExpression}

import scala.concurrent.Future
import scala.util.Try

class BackendWorkflowInitializationActorSpec extends TestKitSuite("BackendWorkflowInitializationActorSpec")
  with FlatSpecLike with Matchers with TableDrivenPropertyChecks {

  behavior of "BackendWorkflowInitializationActorSpec"

  val testPredicateBackendWorkflowInitializationActorRef:
    TestActorRef[TestPredicateBackendWorkflowInitializationActor] =
    TestActorRef[TestPredicateBackendWorkflowInitializationActor]

  val testPredicateBackendWorkflowInitializationActor:
    TestPredicateBackendWorkflowInitializationActor =
    testPredicateBackendWorkflowInitializationActorRef.underlyingActor

  val testContinueOnReturnCode: (Option[WdlValue]) => Boolean = {
    testPredicateBackendWorkflowInitializationActor.continueOnReturnCodePredicate(valueRequired = false)
  }

  val optionalConfig = Option(TestConfig.optionalRuntimeConfig)

  it should "continueOnReturnCodePredicate" in {
    testContinueOnReturnCode(None) should be(true)
    ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(None) should be(true)

    val booleanRows = Table(
      "value",
      true,
      false
    )

    val integerRows = Table(
      "value",
      -1,
      0,
      1,
      1024
    )

    val expressionRows = Table(
      "expression",
      "read_int(\"bad file\")"
    )

    val invalidWdlValueRows = Table(
      "wdlValue",
      WdlString(""),
      WdlString("z"),
      WdlFloat(0.0D),
      WdlArray(WdlArrayType(WdlBooleanType), Seq(WdlBoolean(true))),
      WdlArray(WdlArrayType(WdlFloatType), Seq(WdlFloat(0.0D)))
    )

    forAll(booleanRows) { value =>
      val wdlValue = WdlBoolean(value)
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> wdlValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeFlag(value))
    }

    forAll(booleanRows) { value =>
      val wdlValue = WdlString(value.toString)
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> wdlValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeFlag(value))
    }

    forAll(booleanRows) { value =>
      val wdlValue = WdlExpression.fromString(value.toString)
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      // NOTE: expressions are never valid to validate
    }

    forAll(integerRows) { value =>
      val wdlValue = WdlInteger(value)
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> wdlValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeSet(Set(value)))
    }

    forAll(integerRows) { value =>
      val wdlValue = WdlString(value.toString)
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> wdlValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeSet(Set(value)))
    }

    forAll(integerRows) { value =>
      val wdlValue = WdlExpression.fromString(value.toString)
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      // NOTE: expressions are never valid to validate
    }

    forAll(integerRows) { value =>
      val wdlValue = WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(value)))
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> wdlValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeSet(Set(value)))
    }

    forAll(integerRows) { value =>
      val wdlValue = WdlArray(WdlArrayType(WdlStringType), Seq(WdlString(value.toString)))
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> wdlValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeSet(Set(value)))
    }

    forAll(integerRows) { value =>
      val wdlValue = WdlArray(WdlArrayType(WdlExpressionType), Seq(WdlExpression.fromString(value.toString)))
      val result = false
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      // NOTE: expressions are never valid to validate
    }

    forAll(expressionRows) { expression =>
      val wdlValue = WdlExpression.fromString(expression)
      val result = true
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      // NOTE: expressions are never valid to validate
    }

    forAll(invalidWdlValueRows) { wdlValue =>
      val result = false
      testContinueOnReturnCode(Option(wdlValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalExpression(Option(wdlValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> wdlValue))
      valid.isValid should be(result)
      valid.toEither.left.get.toList should contain theSameElementsAs List(
        "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]"
      )
    }

  }

}

class TestPredicateBackendWorkflowInitializationActor extends BackendWorkflowInitializationActor {

  override val serviceRegistryActor: ActorRef = context.system.deadLetters

  override def calls: Set[TaskCall] = throw new NotImplementedError("calls")

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlValue]) => Boolean] =
    throw new NotImplementedError("runtimeAttributeValidators")

  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WdlValue]] =
    throw new NotImplementedError("coerceDefaultRuntimeAttributes")

  override def beforeAll(): Future[Option[BackendInitializationData]] = throw new NotImplementedError("beforeAll")

  override def validate(): Future[Unit] = throw new NotImplementedError("validate")

  override protected def workflowDescriptor: BackendWorkflowDescriptor =
    throw new NotImplementedError("workflowDescriptor")

  override protected def configurationDescriptor: BackendConfigurationDescriptor = BackendConfigurationDescriptor(TestConfig.sampleBackendRuntimeConfig, ConfigFactory.empty())

  override def continueOnReturnCodePredicate(valueRequired: Boolean)
                                            (wdlExpressionMaybe: Option[WdlValue]): Boolean = {
    super.continueOnReturnCodePredicate(valueRequired)(wdlExpressionMaybe)
  }
}
