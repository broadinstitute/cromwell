package cromwell.backend

import _root_.wdl.draft2.model.types._
import _root_.wdl.draft2.model.WdlExpression
import akka.actor.ActorRef
import akka.testkit.TestActorRef
import com.typesafe.config.ConfigFactory
import cromwell.backend.validation.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet, ContinueOnReturnCodeValidation}
import cromwell.core.{TestKitSuite, WorkflowOptions}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}
import wom.RuntimeAttributesKeys
import wom.expression.WomExpression
import wom.graph.CommandCallNode
import wom.types._
import wom.values._

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

  val testContinueOnReturnCode: (Option[WomValue]) => Boolean = {
    testPredicateBackendWorkflowInitializationActor.continueOnReturnCodePredicate(valueRequired = false)
  }

  val optionalConfig = Option(TestConfig.optionalRuntimeConfig)

  it should "continueOnReturnCodePredicate" in {
    testContinueOnReturnCode(None) should be(true)
    ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(None) should be(true)

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
      "womValue",
      WomString(""),
      WomString("z"),
      WomFloat(0.0D),
      WomArray(WomArrayType(WomBooleanType), Seq(WomBoolean(true))),
      WomArray(WomArrayType(WomFloatType), Seq(WomFloat(0.0D)))
    )

    forAll(booleanRows) { value =>
      val womValue = WomBoolean(value)
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> womValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeFlag(value))
    }

    forAll(booleanRows) { value =>
      val womValue = WomString(value.toString)
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> womValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeFlag(value))
    }

    forAll(booleanRows) { value =>
      val womValue = WdlExpression.fromString(value.toString)
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      // NOTE: expressions are never valid to validate
    }

    forAll(integerRows) { value =>
      val womValue = WomInteger(value)
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> womValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeSet(Set(value)))
    }

    forAll(integerRows) { value =>
      val womValue = WomString(value.toString)
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> womValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeSet(Set(value)))
    }

    forAll(integerRows) { value =>
      val womValue = WdlExpression.fromString(value.toString)
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      // NOTE: expressions are never valid to validate
    }

    forAll(integerRows) { value =>
      val womValue = WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(value)))
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> womValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeSet(Set(value)))
    }

    forAll(integerRows) { value =>
      val womValue = WomArray(WomArrayType(WomStringType), Seq(WomString(value.toString)))
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> womValue))
      valid.isValid should be(result)
      valid.toEither.right.get should be(ContinueOnReturnCodeSet(Set(value)))
    }

    forAll(integerRows) { value =>
      val womValue = WomArray(WomArrayType(WdlExpressionType), Seq(WdlExpression.fromString(value.toString)))
      val result = false
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      // NOTE: expressions are never valid to validate
    }

    forAll(expressionRows) { expression =>
      val womValue = WdlExpression.fromString(expression)
      val result = true
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      // NOTE: expressions are never valid to validate
    }

    forAll(invalidWdlValueRows) { womValue =>
      val result = false
      testContinueOnReturnCode(Option(womValue)) should be(result)
      ContinueOnReturnCodeValidation.default(optionalConfig).validateOptionalWomValue(Option(womValue)) should be(result)
      val valid =
        ContinueOnReturnCodeValidation.default(optionalConfig).validate(Map(RuntimeAttributesKeys.ContinueOnReturnCodeKey -> womValue))
      valid.isValid should be(result)
      valid.toEither.left.get.toList should contain theSameElementsAs List(
        "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]"
      )
    }

  }

}

class TestPredicateBackendWorkflowInitializationActor extends BackendWorkflowInitializationActor {

  override val serviceRegistryActor: ActorRef = context.system.deadLetters

  override def calls: Set[CommandCallNode] = throw new UnsupportedOperationException("calls")

  override protected def runtimeAttributeValidators: Map[String, (Option[WomExpression]) => Boolean] =
    throw new UnsupportedOperationException("runtimeAttributeValidators")

  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WomValue]] =
    throw new UnsupportedOperationException("coerceDefaultRuntimeAttributes")

  override def beforeAll(): Future[Option[BackendInitializationData]] = throw new UnsupportedOperationException("beforeAll")

  override def validate(): Future[Unit] = throw new UnsupportedOperationException("validate")

  override protected def workflowDescriptor: BackendWorkflowDescriptor =
    throw new UnsupportedOperationException("workflowDescriptor")

  override protected def configurationDescriptor: BackendConfigurationDescriptor = BackendConfigurationDescriptor(TestConfig.sampleBackendRuntimeConfig, ConfigFactory.empty())

  override def continueOnReturnCodePredicate(valueRequired: Boolean)
                                            (wdlExpressionMaybe: Option[WomValue]): Boolean = {
    super.continueOnReturnCodePredicate(valueRequired)(wdlExpressionMaybe)
  }
}
