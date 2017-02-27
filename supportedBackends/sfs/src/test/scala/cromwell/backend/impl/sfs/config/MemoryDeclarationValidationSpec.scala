package cromwell.backend.impl.sfs.config

import com.typesafe.config.ConfigFactory
import cromwell.backend.MemorySize
import cromwell.backend.validation.{RuntimeAttributesKeys, ValidatedRuntimeAttributes}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.parser.MemoryUnit
import wdl4s.values.{WdlFloat, WdlInteger}

class MemoryDeclarationValidationSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  behavior of "MemoryDeclarationValidation"

  val validDeclaredAmounts = Table(
    ("declaration", "runtimeAmount", "expectedDefaultAmount", "expectedExtracted"),
    ("Int memory", Option(2), None, Option(WdlInteger(2 * 1000 * 1000 * 1000))),
    ("Int memory_gb", Option(2), None, Option(WdlInteger(2))),
    ("Int memory_gb = 3", None, Option(3), None),
    ("Int memory_gb = 3", Option(2), Option(3), Option(WdlInteger(2))),
    ("Int? memory_gb", None, None, None),
    ("Int? memory_gb", Option(2), None, Option(WdlInteger(2))),
    ("Int? memory_gb = 3", None, Option(3), None),
    ("Int? memory_gb = 3", Option(2), Option(3), Option(WdlInteger(2))),
    ("Float memory", Option(2), None, Option(WdlFloat(2 * 1000 * 1000 * 1000))),
    ("Float memory_gb", Option(2), None, Option(WdlFloat(2))),
    ("Float memory_gb = 3.0", None, Option(3), None),
    ("Float memory_gb = 3.0", Option(2), Option(3), Option(WdlFloat(2))),
    ("Float? memory_gb", None, None, None),
    ("Float? memory_gb", Option(2), None, Option(WdlFloat(2))),
    ("Float? memory_gb = 3.0", None, Option(3), None),
    ("Float? memory_gb = 3.0", Option(2), Option(3), Option(WdlFloat(2)))
  )

  forAll(validDeclaredAmounts) { (declaration, runtimeAmount, expectedDefaultAmount, expectedExtracted) =>
    it should s"extract memory from declared $declaration with memory set to ${runtimeAmount.getOrElse("none")}" in {
      val config = ConfigFactory.parseString(
        s"""|submit = "anything"
            |${ConfigConstants.RuntimeAttributesConfig} = "$declaration"
            |""".stripMargin)

      val configWdlNamespace = new ConfigWdlNamespace(config)
      val runtimeDeclaration = configWdlNamespace.runtimeDeclarations.head
      val memoryDeclarationValidation = new MemoryDeclarationValidation(runtimeDeclaration)
      val attributes = runtimeAmount
        .map(amount => RuntimeAttributesKeys.MemoryKey -> MemorySize(amount.toDouble, MemoryUnit.GB))
        .toMap
      val validatedRuntimeAttributes = ValidatedRuntimeAttributes(attributes)

      val default = memoryDeclarationValidation.makeValidation().runtimeAttributeDefinition.factoryDefault
      val extracted = memoryDeclarationValidation.extractWdlValueOption(validatedRuntimeAttributes)

      val expectedDefault = expectedDefaultAmount
        .map(amount => WdlInteger(MemorySize(amount.toDouble, MemoryUnit.GB).bytes.toInt))

      MemoryDeclarationValidation.isMemoryDeclaration(runtimeDeclaration.unqualifiedName) should be(true)
      default should be(expectedDefault)
      extracted should be(expectedExtracted)
    }
  }

  val badSyntaxDeclarations = Table(
    "declaration",
    "Int memory_gb = 3.0",
    "Float memory_gb = 3"
  )

  forAll(badSyntaxDeclarations) { declaration =>
    it should s"throw a syntax error for memory declaration $declaration" in {
      val config = ConfigFactory.parseString(
        s"""|submit = "anything"
            |${ConfigConstants.RuntimeAttributesConfig} = "$declaration"
            |""".stripMargin)

      val expectedException = intercept[RuntimeException](new ConfigWdlNamespace(config))
      expectedException.getMessage should startWith("Error parsing generated wdl:\n")
    }
  }

  val invalidDeclarations = Table(
    "declaration",
    "Int mem",
    "Int memory_badunit",
    "Float memory_badunit"
  )

  forAll(invalidDeclarations) { declaration =>
    it should s"not identify $declaration as a memory declaration" in {
      val config = ConfigFactory.parseString(
        s"""|submit = "anything"
            |${ConfigConstants.RuntimeAttributesConfig} = "$declaration"
            |""".stripMargin)

      val configWdlNamespace = new ConfigWdlNamespace(config)
      val runtimeDeclaration = configWdlNamespace.runtimeDeclarations.head
      MemoryDeclarationValidation.isMemoryDeclaration(runtimeDeclaration.unqualifiedName) should be(false)
    }
  }
}
