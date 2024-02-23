package cromwell.backend.impl.sfs.config

import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.impl.sfs.config.ConfigConstants.{MemoryRuntimeAttribute, _}
import cromwell.backend.validation.ValidatedRuntimeAttributes
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.values.{WomFloat, WomLong}

class MemoryDeclarationValidationSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with TableDrivenPropertyChecks {
  behavior of "MemoryDeclarationValidation"

  val validDeclaredAmounts = Table(
    ("declaration", "runtimeAmount", "expectedDefaultAmount", "expectedExtracted"),
    ("Int memory", Option(2), None, Option(WomLong(2L * 1024L * 1024L * 1024L))),
    ("Int memory_gb", Option(2), None, Option(WomLong(2))),
    ("Int memory_gb = 3", None, Option(3), None),
    ("Int memory_gb = 3", Option(2), Option(3), Option(WomLong(2))),
    ("Int? memory_gb", None, None, None),
    ("Int? memory_gb", Option(2), None, Option(WomLong(2))),
    ("Int? memory_gb = 3", None, Option(3), None),
    ("Int? memory_gb = 3", Option(2), Option(3), Option(WomLong(2))),
    ("Float memory", Option(2), None, Option(WomFloat(2f * 1024f * 1024f * 1024f))),
    ("Float memory_gb", Option(2), None, Option(WomFloat(2))),
    ("Float memory_gb = 3.0", None, Option(3), None),
    ("Float memory_gb = 3.0", Option(2), Option(3), Option(WomFloat(2))),
    ("Float memory_gb = 3", None, Option(3), None),
    ("Float memory_gb = 3", Option(2), Option(3), Option(WomFloat(2))),
    ("Float? memory_gb", None, None, None),
    ("Float? memory_gb", Option(2), None, Option(WomFloat(2))),
    ("Float? memory_gb = 3.0", None, Option(3), None),
    ("Float? memory_gb = 3.0", Option(2), Option(3), Option(WomFloat(2))),
    ("Float? memory_gb = 3", None, Option(3), None),
    ("Float? memory_gb = 3", Option(2), Option(3), Option(WomFloat(2)))
  )

  forAll(validDeclaredAmounts) { (declaration, runtimeAmount, expectedDefaultAmount, expectedExtracted) =>
    it should s"extract memory from declared $declaration with memory set to ${runtimeAmount.getOrElse("none")}" in {
      val memoryPrefix =
        if (declaration.contains("Min")) MemoryMinRuntimeAttributePrefix
        else if (declaration.contains("Max")) MemoryMaxRuntimeAttributePrefix
        else MemoryRuntimeAttributePrefix

      val config = ConfigFactory.parseString(s"""|submit = "anything"
                                                 |${ConfigConstants.RuntimeAttributesConfig} = "$declaration"
                                                 |""".stripMargin)

      val configWdlNamespace = new ConfigWdlNamespace(config)
      val runtimeDeclaration = configWdlNamespace.runtimeDeclarations.head
      val memoryDeclarationValidation =
        new MemoryDeclarationValidation(runtimeDeclaration, MemoryRuntimeAttribute, memoryPrefix)
      val attributes = runtimeAmount
        .map(amount => MemoryRuntimeAttribute -> MemorySize(amount.toDouble, MemoryUnit.GB))
        .toMap
      val validatedRuntimeAttributes = ValidatedRuntimeAttributes(attributes)

      val default = memoryDeclarationValidation.makeValidation().runtimeAttributeDefinition.factoryDefault
      val extracted = memoryDeclarationValidation.extractWdlValueOption(validatedRuntimeAttributes)

      val expectedDefault = expectedDefaultAmount
        .map(amount => WomLong(MemorySize(amount.toDouble, MemoryUnit.GB).bytes.toLong))

      MemoryDeclarationValidation.isMemoryDeclaration(runtimeDeclaration.unqualifiedName,
                                                      MemoryRuntimeAttribute,
                                                      memoryPrefix
      ) should be(true)
      default should be(expectedDefault)
      extracted should be(expectedExtracted)
    }
  }

  val badSyntaxDeclarations = Table(
    "declaration",
    "Int memory_gb = 3.0"
  )

  forAll(badSyntaxDeclarations) { declaration =>
    it should s"throw a syntax error for memory declaration $declaration" in {
      val config = ConfigFactory.parseString(s"""|submit = "anything"
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
      val config = ConfigFactory.parseString(s"""|submit = "anything"
                                                 |${ConfigConstants.RuntimeAttributesConfig} = "$declaration"
                                                 |""".stripMargin)

      val configWdlNamespace = new ConfigWdlNamespace(config)
      val runtimeDeclaration = configWdlNamespace.runtimeDeclarations.head
      MemoryDeclarationValidation.isMemoryDeclaration(runtimeDeclaration.unqualifiedName,
                                                      MemoryRuntimeAttribute,
                                                      MemoryRuntimeAttributePrefix
      ) should be(false)
    }
  }
}
