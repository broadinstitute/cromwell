package cromwell.backend.impl.sfs.config

import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.{TableDrivenPropertyChecks, TableFor4}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.values.{WomString, WomValue}
import cats.syntax.validated._
import common.validation.ErrorOr._

class ConfigInitializationActorSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with TableDrivenPropertyChecks {

  it should "read the list of call cache attributes from config" in {
    val tripleQuote = "\"\"\""
    val singularityBackendConfigString =
      s"""run-in-background = true
         | runtime-attributes = ${tripleQuote}
         |       String singularity
         |       Int uncacheworthy1
         |       Int uncacheworthy2
         |     ${tripleQuote}
         | runtime-attributes-for-caching = {
         |   "singularity": true
         |   "uncacheworthy1": false
         |   # NB: no specific entry for uncacheworthy2 - the default should already be 'false'.
         | }
         | submit-docker = ${tripleQuote}
         |       singularity exec --bind $${cwd}:$${docker_cwd} docker://$${singularity} $${job_shell} $${script}
         |     ${tripleQuote}
         |""".stripMargin

    val singularityBackendConfig = ConfigFactory.parseString(singularityBackendConfigString)

    // Mirroring how the declarations are made in ConfigInitializationActor:
    lazy val configWdlNamespace = new ConfigWdlNamespace(singularityBackendConfig)
    lazy val declarationValidations: Seq[DeclarationValidation] =
      DeclarationValidation.fromDeclarations(configWdlNamespace.runtimeDeclarations,
                                             configWdlNamespace.callCachedRuntimeAttributes
      )

    declarationValidations.exists(p => p.key == "singularity" && p.makeValidation().usedInCallCaching) should be(true)
    declarationValidations.exists(p => p.key == "uncacheworthy1" && p.makeValidation().usedInCallCaching) should be(
      false
    )
    declarationValidations.exists(p => p.key == "uncacheworthy2" && p.makeValidation().usedInCallCaching) should be(
      false
    )

  }

  private val attributeValidationTests: TableFor4[String, String, Map[String, WomValue], ErrorOr[_]] = Table(
    ("description", "declaration", "womValue", "expected"),
    ("present required custom", "Int my_attr", Map("my_attr" -> WomString("1")), 1.valid),
    ("missing required custom",
     "Int my_attr",
     Map(),
     "Expecting my_attr runtime attribute to be an Integer".invalidNel
    ),
    ("present defaulted custom", "Int my_attr = 415", Map("my_attr" -> WomString("1")), 1.valid),
    ("missing defaulted custom", "Int my_attr = 415", Map(), 415.valid),
    ("present optional custom", "Int? my_attr", Map("my_attr" -> WomString("1")), Option(1).valid),
    ("missing optional custom", "Int? my_attr", Map(), None.valid),
    ("present defaulted optional custom", "Int? my_attr = 415", Map("my_attr" -> WomString("1")), Option(1).valid),
    ("missing defaulted optional custom", "Int? my_attr = 415", Map(), Option(415).valid),
    ("present required memory",
     "Float memory_gb",
     Map("memory" -> WomString("3 GB")),
     MemorySize(3, MemoryUnit.GB).valid
    ),
    ("missing required memory",
     "Float memory_gb",
     Map(),
     "Expecting memory runtime attribute to be an Integer or String with format '8 GB'. Exception: Not supported WDL type value".invalidNel
    ),
    ("present defaulted memory",
     "Float memory_gb = 1",
     Map("memory" -> WomString("3 GB")),
     MemorySize(3, MemoryUnit.GB).valid
    ),
    ("missing defaulted memory", "Float memory_gb = 1", Map(), MemorySize(1, MemoryUnit.GB).valid),
    ("present optional memory",
     "Float? memory_gb",
     Map("memory" -> WomString("3 GB")),
     Option(MemorySize(3, MemoryUnit.GB)).valid
    ),
    ("missing optional memory", "Float? memory_gb", Map(), None.valid),
    ("present defaulted optional memory",
     "Float? memory_gb = 1",
     Map("memory" -> WomString("3 GB")),
     Option(MemorySize(3, MemoryUnit.GB)).valid
    ),
    ("missing defaulted optional memory", "Float? memory_gb = 1", Map(), Option(MemorySize(1, MemoryUnit.GB)).valid)
  )
  forAll(attributeValidationTests) { (description, declaration, values, expected) =>
    it should s"validate a $description attribute" in {
      val backendConfig = ConfigFactory.parseString(
        s"""submit: "unused"
           |runtime-attributes: "$declaration"
           |""".stripMargin
      )
      lazy val configWdlNamespace = new ConfigWdlNamespace(backendConfig)
      val validations = DeclarationValidation.fromDeclarations(
        configWdlNamespace.runtimeDeclarations,
        configWdlNamespace.callCachedRuntimeAttributes
      )
      val validation = validations.head.makeValidation()
      val actual = validation.validate(values)
      actual should be(expected)
    }
  }
}
