package languages.wdl.draft3

import com.typesafe.config.{Config, ConfigFactory}
import common.validation.ErrorOr.ErrorOr
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesWithoutImports}
import cromwell.languages.LanguageFactory
import cromwell.languages.util.ImportResolver
import languages.wdl.draft3.WdlDraft3CachingSpec.EvaluationCountingDraft3Factory
import org.scalatest.{FlatSpec, Matchers}
import wom.ResolvedImportRecord
import wom.core.{WorkflowOptionsJson, WorkflowSource}
import wom.executable.WomBundle
import wom.expression.NoIoFunctionSet

class WdlDraft3CachingSpec extends FlatSpec with Matchers {

  val languageConfig = ConfigFactory.parseString(
    """{
      |  strict-validation: true
      |  enabled: true
      |  caching {
      |    enabled: true
      |    ttl: 3 minutes
      |    size: 50
      |    concurrency: 9
      |  }
      |}
      |""".stripMargin

  )



  it should "only evaluate files once" in {
    val invalidWorkflowSource =
      """
        |blah blah this isn't a real workflow
      """.stripMargin

    val validWorkflowSource =
      """version 1.0
        |
        |task hello {
        |  input {
        |    String addressee
        |  }
        |  command {
        |    echo "Hello ${addressee}!"
        |  }
        |  runtime {
        |      docker: "ubuntu:latest"
        |  }
        |  output {
        |    String salutation = read_string(stdout())
        |  }
        |}
        |
        |workflow wf_hello {
        |  String wf_hello_input = "world"
        |
        |  call hello { input: addressee = wf_hello_input }
        |
        |  output {
        |    String salutation = hello.salutation
        |  }
        |}
      """.stripMargin

    val factory = new EvaluationCountingDraft3Factory(languageConfig)

    def validate(workflowSource: WorkflowSource) = factory.validateNamespace(
      WorkflowSourceFilesWithoutImports(
        Option(workflowSource),
        None,
        None,
        None,
        None,
        "{}",
        WorkflowOptions.empty,
        "{}",
        workflowOnHold = false,
        Seq.empty
      ),
      workflowSource,
      WorkflowOptions.empty,
      importLocalFilesystem = false,
      WorkflowId.randomId(),
      NoIoFunctionSet,
      List.empty
    )

    // Check the valid workflow twice:
    validate(validWorkflowSource).isRight.unsafeRunSync() should be(true)
    validate(validWorkflowSource).isRight.unsafeRunSync() should be(true)

    // But we only evaludated it once:
    factory.evaluationCount should be(1)

    // Check the invalid workflow twice:
    validate(invalidWorkflowSource).isRight.unsafeRunSync() should be(false)
    validate(invalidWorkflowSource).isRight.unsafeRunSync() should be(false)

    // The factory only ran one extra evaluation:
    factory.evaluationCount should be(2)

    // Run over the two workflows a few more times:
    validate(validWorkflowSource).isRight.unsafeRunSync() should be(true)
    validate(validWorkflowSource).isRight.unsafeRunSync() should be(true)
    validate(invalidWorkflowSource).isRight.unsafeRunSync() should be(false)
    validate(invalidWorkflowSource).isRight.unsafeRunSync() should be(false)

    // No additional evaluations were needed:
    factory.evaluationCount should be(2)
  }

}

object WdlDraft3CachingSpec {
  class EvaluationCountingDraft3Factory(languageConfig: Config) extends WdlDraft3LanguageFactory(languageConfig) {

    var evaluationCount = 0

    override protected def makeWomBundle(workflowSource: WorkflowSource, workflowSourceOrigin: Option[ResolvedImportRecord], workflowOptionsJson: WorkflowOptionsJson, importResolvers: List[ImportResolver.ImportResolver], languageFactories: List[LanguageFactory], convertNestedScatterToSubworkflow: Boolean): ErrorOr[WomBundle] = {
      evaluationCount = evaluationCount + 1
      super.makeWomBundle(workflowSource, workflowSourceOrigin, workflowOptionsJson, importResolvers, languageFactories, convertNestedScatterToSubworkflow)
    }

  }
}
