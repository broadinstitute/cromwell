package languages.wdl.draft2

import cats.instances.list._
import cats.syntax.functor._
import com.typesafe.config.ConfigFactory
import common.Checked
import cromwell.core.{CacheConfig, WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.languages.util.ImportResolver.HttpResolver
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import wom.expression.NoIoFunctionSet

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps

class NamespaceCacheSpec extends FlatSpec with BeforeAndAfterAll with Matchers {

  "The WdlDraft2LanguageFactory" should "parse a fully enabled config correctly" in {
    val factory = new WdlDraft2LanguageFactory(ConfigFactory.parseString(EnabledConfig))
    factory.cacheConfig.shouldBe(Option(CacheConfig(ttl = 3 minutes, size = 50L, concurrency = 9)))
  }

  it should "not make a cache config for a disabled language factory config" in {
    val factory = new WdlDraft2LanguageFactory(ConfigFactory.parseString(DisabledLanguageConfig))
    factory.cacheConfig.shouldBe(None)
  }

  it should "not make a cache config for an enabled language factory config with disabled caching" in {
    val factory = new WdlDraft2LanguageFactory(ConfigFactory.parseString(DisabledCacheConfig))
    factory.cacheConfig.shouldBe(None)
  }

  it should "use a fully enabled namespace cache" in {
    validate(EnabledConfig, (1 to 10).toList.as(1))
  }

  it should "not use a disabled namespace cache" in {
    validate(DisabledCacheConfig, (1 to 3).toList)
  }

  def validate(configString: String, expectations: List[Int]): Unit = {
    val config = ConfigFactory.parseString(configString)
    val factory = new WdlDraft2LanguageFactory(config)
    val EmptyJson = "{}"
    val collection = WorkflowSourceFilesCollection(
      workflowSource = Option(ThreeStep),
      workflowUrl = None,
      workflowRoot = None,
      workflowType = None,
      workflowTypeVersion = None,
      inputsJson = EmptyJson,
      workflowOptions = WorkflowOptions.empty,
      labelsJson = EmptyJson,
      importsFile = None,
      workflowOnHold = false,
      warnings = Seq.empty
    )

    var lookupCount = 0
    val countingResolver = new HttpResolver() {
      override def pathToLookup(str: String): Checked[String] = {
        lookupCount = lookupCount + 1
        super.pathToLookup(str)
      }
    }

    def validate = {
      val futureNamespace = factory.validateNamespace(
        source = collection,
        workflowSource = ThreeStep,
        workflowOptions = WorkflowOptions(new spray.json.JsObject(Map.empty)),
        importLocalFilesystem = false,
        workflowIdForLogging = WorkflowId.randomId(),
        ioFunctions = NoIoFunctionSet,
        importResolvers = List(countingResolver)).value.unsafeToFuture()
      Await.result(futureNamespace, Duration.Inf).right.get
    }

    expectations foreach { e =>
      validate
      lookupCount shouldBe e
    }
  }

  val ThreeStep =
    """
      |# Not actually used by the workflow, just here to exercise the Importer.
      |import "https://raw.githubusercontent.com/broadinstitute/cromwell/f2d1c3bdd5535d9f6a997eadaf136742d86adbe5/centaur/src/main/resources/standardTestCases/aliased_subworkflows/subworkflow.wdl"
      |
      |task hello {
      |  command {
      |    echo "Hello esteemed user!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |}
      |
      |workflow wf_hello {
      |  call hello
      |}""".stripMargin

  val EnabledConfig =
    """
      |{
      |  enabled: true
      |  strict-validation: true
      |  caching {
      |    enabled: true
      |    ttl: 3 minutes
      |    size: 50
      |    concurrency: 9
      |  }
      |}
    """.stripMargin

  val DisabledCacheConfig =
    """
      |{
      |  enabled: true
      |  strict-validation: true
      |  caching {
      |    enabled: false
      |    ttl: 3 minutes
      |    size: 50
      |    concurrency: 9
      |  }
      |}
    """.stripMargin

  val DisabledLanguageConfig =
    """
      |{
      |  enabled: false
      |  strict-validation: true
      |  caching {
      |    enabled: true
      |    ttl: 3 minutes
      |    size: 50
      |    concurrency: 9
      |  }
      |}
    """.stripMargin
}
