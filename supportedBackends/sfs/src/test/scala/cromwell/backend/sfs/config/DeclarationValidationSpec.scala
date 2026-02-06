package cromwell.backend.sfs.config

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.impl.sfs.config.DeclarationValidation
import cromwell.backend.validation.{Containers, ValidatedRuntimeAttributes}
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import wdl.draft2.model._
import wom.types.{WomIntegerType, WomOptionalType, WomStringType}
import wom.values.{WomInteger, WomString}

class DeclarationValidationSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with TableDrivenPropertyChecks {
  behavior of "DeclarationValidation"

  def validateCpu(key: String) = {
    val expression = WdlExpression.fromString("5")
    val declarationValidation = DeclarationValidation.fromDeclaration(callCachedRuntimeAttributesMap = Map.empty)(
      Declaration(WomIntegerType, key, Option(expression), None, null)
    )
    declarationValidation.extractWdlValueOption(
      ValidatedRuntimeAttributes(Map(key -> refineMV[Positive](5)))
    ) shouldBe Some(WomInteger(5))
  }

  it should "validate cpu attributes" in {
    val keys = Table(
      "key",
      "cpu"
    )

    forAll(keys)(validateCpu)
  }

  def validateContainer(key: String) = {
    val expression = WdlExpression.fromString("\"ubuntu:latest\"")
    val declarationValidation = DeclarationValidation.fromDeclaration(callCachedRuntimeAttributesMap = Map.empty)(
      Declaration(WomStringType, key, Option(expression), None, null)
    )
    declarationValidation.extractWdlValueOption(
      ValidatedRuntimeAttributes(Map(key -> Containers("ubuntu:latest")))
    ) shouldBe Some(WomString("ubuntu:latest"))
  }

  it should "validate docker attribute with Containers type" in {
    validateContainer("docker")
  }

  it should "validate container attribute with Containers type" in {
    validateContainer("container")
  }

  it should "return None for absent optional docker attribute" in {
    val expression = WdlExpression.fromString("\"ubuntu:latest\"")
    val declarationValidation = DeclarationValidation.fromDeclaration(callCachedRuntimeAttributesMap = Map.empty)(
      Declaration(WomOptionalType(WomStringType), "docker", Option(expression), None, null)
    )
    declarationValidation.extractWdlValueOption(
      ValidatedRuntimeAttributes(Map.empty)
    ) shouldBe None
  }

  it should "extract first image when multiple container images are provided" in {
    val expression = WdlExpression.fromString("\"ubuntu:latest\"")
    val declarationValidation = DeclarationValidation.fromDeclaration(callCachedRuntimeAttributesMap = Map.empty)(
      Declaration(WomStringType, "docker", Option(expression), None, null)
    )
    declarationValidation.extractWdlValueOption(
      ValidatedRuntimeAttributes(Map("docker" -> Containers(List("first:1.0", "second:2.0", "third:3.0"))))
    ) shouldBe Some(WomString("first:1.0"))
  }
}
