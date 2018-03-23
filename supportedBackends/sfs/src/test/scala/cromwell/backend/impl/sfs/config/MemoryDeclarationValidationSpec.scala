package cromwell.backend.impl.sfs.config

import com.typesafe.config.ConfigFactory
import cromwell.backend.impl.sfs.config.ConfigConstants.{MemoryRuntimeAttribute, _}
import cromwell.backend.validation.ValidatedRuntimeAttributes
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.values.{WomFloat, WomInteger}

class MemoryDeclarationValidationSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  behavior of "MemoryDeclarationValidation"

  val validDeclaredAmounts = Table(
    ("declaration", "runtimeAmount", "expectedDefaultAmount", "expectedExtracted"),
    ("Int memory", Option(2), None, Option(WomInteger(2 * 1000 * 1000 * 1000))),
    ("Int memory_gb", Option(2), None, Option(WomInteger(2))),
    ("Int memory_gb = 3", None, Option(3), None),
    ("Int memory_gb = 3", Option(2), Option(3), Option(WomInteger(2))),
    ("Int? memory_gb", None, None, None),
    ("Int? memory_gb", Option(2), None, Option(WomInteger(2))),
    ("Int? memory_gb = 3", None, Option(3), None),
    ("Int? memory_gb = 3", Option(2), Option(3), Option(WomInteger(2))),
    ("Float memory", Option(2), None, Option(WomFloat(2 * 1000 * 1000 * 1000))),
    ("Float memory_gb", Option(2), None, Option(WomFloat(2))),
    ("Float memory_gb = 3.0", None, Option(3), None),
    ("Float memory_gb = 3.0", Option(2), Option(3), Option(WomFloat(2))),
    ("Float? memory_gb", None, None, None),
    ("Float? memory_gb", Option(2), None, Option(WomFloat(2))),
    ("Float? memory_gb = 3.0", None, Option(3), None),
    ("Float? memory_gb = 3.0", Option(2), Option(3), Option(WomFloat(2))),
    ("Int memoryMin", Option(2), None, Option(WomInteger(2 * 1000 * 1000 * 1000))),
    ("Int memoryMin_gb", Option(2), None, Option(WomInteger(2))),
    ("Int memoryMin_gb = 3", None, Option(3), None),
    ("Int memoryMin_gb = 3", Option(2), Option(3), Option(WomInteger(2))),
    ("Int? memoryMin_gb", None, None, None),
    ("Int? memoryMin_gb", Option(2), None, Option(WomInteger(2))),
    ("Int? memoryMin_gb = 3", None, Option(3), None),
    ("Int? memoryMin_gb = 3", Option(2), Option(3), Option(WomInteger(2))),
    ("Float memoryMin", Option(2), None, Option(WomFloat(2 * 1000 * 1000 * 1000))),
    ("Float memoryMin_gb", Option(2), None, Option(WomFloat(2))),
    ("Float memoryMin_gb = 3.0", None, Option(3), None),
    ("Float memoryMin_gb = 3.0", Option(2), Option(3), Option(WomFloat(2))),
    ("Float? memoryMin_gb", None, None, None),
    ("Float? memoryMin_gb", Option(2), None, Option(WomFloat(2))),
    ("Float? memoryMin_gb = 3.0", None, Option(3), None),
    ("Float? memoryMin_gb = 3.0", Option(2), Option(3), Option(WomFloat(2))),
    ("Int memoryMax", Option(2), None, Option(WomInteger(2 * 1000 * 1000 * 1000))),
    ("Int memoryMax_gb", Option(2), None, Option(WomInteger(2))),
    ("Int memoryMax_gb = 3", None, Option(3), None),
    ("Int memoryMax_gb = 3", Option(2), Option(3), Option(WomInteger(2))),
    ("Int? memoryMax_gb", None, None, None),
    ("Int? memoryMax_gb", Option(2), None, Option(WomInteger(2))),
    ("Int? memoryMax_gb = 3", None, Option(3), None),
    ("Int? memoryMax_gb = 3", Option(2), Option(3), Option(WomInteger(2))),
    ("Float memoryMax", Option(2), None, Option(WomFloat(2 * 1000 * 1000 * 1000))),
    ("Float memoryMax_gb", Option(2), None, Option(WomFloat(2))),
    ("Float memoryMax_gb = 3.0", None, Option(3), None),
    ("Float memoryMax_gb = 3.0", Option(2), Option(3), Option(WomFloat(2))),
    ("Float? memoryMax_gb", None, None, None),
    ("Float? memoryMax_gb", Option(2), None, Option(WomFloat(2))),
    ("Float? memoryMax_gb = 3.0", None, Option(3), None),
    ("Float? memoryMax_gb = 3.0", Option(2), Option(3), Option(WomFloat(2)))
  )

  forAll(validDeclaredAmounts) { (declaration, runtimeAmount, expectedDefaultAmount, expectedExtracted) =>
    it should s"extract memory from declared $declaration with memory set to ${runtimeAmount.getOrElse("none")}" in {
      val memoryKey = if (declaration.contains("Min")) MemoryMinRuntimeAttribute
      else if (declaration.contains("Max")) MemoryMaxRuntimeAttribute
      else MemoryRuntimeAttribute

      val memoryPrefix = if (declaration.contains("Min")) MemoryMinRuntimeAttributePrefix
      else if (declaration.contains("Max")) MemoryMaxRuntimeAttributePrefix
      else MemoryRuntimeAttributePrefix
      
      val config = ConfigFactory.parseString(
        s"""|submit = "anything"
            |${ConfigConstants.RuntimeAttributesConfig} = "$declaration"
            |""".stripMargin)

      val configWdlNamespace = new ConfigWdlNamespace(config)
      val runtimeDeclaration = configWdlNamespace.runtimeDeclarations.head
      val memoryDeclarationValidation = new MemoryDeclarationValidation(runtimeDeclaration,
        memoryKey, memoryPrefix)
      val attributes = runtimeAmount
        .map(amount => memoryKey -> MemorySize(amount.toDouble, MemoryUnit.GB))
        .toMap
      val validatedRuntimeAttributes = ValidatedRuntimeAttributes(attributes)

      val default = memoryDeclarationValidation.makeValidation().runtimeAttributeDefinition.factoryDefault
      val extracted = memoryDeclarationValidation.extractWdlValueOption(validatedRuntimeAttributes)

      val expectedDefault = expectedDefaultAmount
        .map(amount => WomInteger(MemorySize(amount.toDouble, MemoryUnit.GB).bytes.toInt))

      MemoryDeclarationValidation.isMemoryDeclaration(runtimeDeclaration.unqualifiedName,
        memoryKey, memoryPrefix) should be(true)
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
      MemoryDeclarationValidation.isMemoryDeclaration(runtimeDeclaration.unqualifiedName,
        MemoryRuntimeAttribute, MemoryRuntimeAttributePrefix) should be(false)
    }
  }
}
