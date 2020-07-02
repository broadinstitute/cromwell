package cromwell.backend.impl.sfs.config

import com.typesafe.config.ConfigFactory
import org.scalatest.{FlatSpec, Matchers}

class ConfigInitializationActorSpec extends FlatSpec with Matchers {

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
    lazy val declarationValidations: Seq[DeclarationValidation] = {
      DeclarationValidation.fromDeclarations(configWdlNamespace.runtimeDeclarations, configWdlNamespace.callCachedRuntimeAttributes)
    }

    declarationValidations.exists(p => p.key == "singularity" && p.makeValidation().usedInCallCaching) should be(true)
    declarationValidations.exists(p => p.key == "uncacheworthy1" && p.makeValidation().usedInCallCaching) should be(false)
    declarationValidations.exists(p => p.key == "uncacheworthy2" && p.makeValidation().usedInCallCaching) should be(false)

  }

}

