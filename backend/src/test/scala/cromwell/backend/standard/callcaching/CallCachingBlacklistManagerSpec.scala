package cromwell.backend.standard.callcaching

import akka.event.NoLogging
import com.typesafe.config.ConfigFactory
import cromwell.core._
import org.scalatest.{FlatSpec, Matchers}
import spray.json._


class CallCachingBlacklistManagerSpec extends FlatSpec with Matchers {
  behavior of "CallCachingBlacklistManager"

  //noinspection RedundantDefaultArgument
  val workflowSourcesNoGrouping = WorkflowSourceFilesWithoutImports(
    workflowSource = None,
    workflowUrl = None,
    workflowRoot = None,
    workflowType = None,
    workflowTypeVersion = None,
    inputsJson = "",
    workflowOptions = WorkflowOptions(JsObject.empty),
    labelsJson = "",
    workflowOnHold = false,
    warnings = List.empty
  )

  val workflowSourcesYesGrouping = workflowSourcesNoGrouping.copy(
    workflowOptions = WorkflowOptions(""" { "google_project": "blacklist_group_testing" } """.parseJson.asJsObject)
  )

  val workflowNoGrouping = new HasWorkflowIdAndSources {
    override def sources: WorkflowSourceFilesCollection = workflowSourcesNoGrouping
    override def id: WorkflowId = WorkflowId.randomId()
  }

  val workflowYesGrouping1 = new HasWorkflowIdAndSources {
    override def sources: WorkflowSourceFilesCollection = workflowSourcesYesGrouping
    override def id: WorkflowId = WorkflowId.randomId()
  }

  val workflowYesGrouping2 = new HasWorkflowIdAndSources {
    override def sources: WorkflowSourceFilesCollection = workflowSourcesYesGrouping
    override def id: WorkflowId = WorkflowId.randomId()
  }

  it should "be off by default" in {
    val configString = ""
    val manager = new CallCachingBlacklistManager(ConfigFactory.parseString(configString), logger = NoLogging)

    manager.blacklistCacheFor(workflowNoGrouping) shouldBe None
  }

  it should "be on with default values if blacklisting is enabled" in {
    val configString =
      """
        |call-caching {
        |  blacklist-cache {
        |    enabled: true
        |  }
        |}
        |""".stripMargin
    val manager = new CallCachingBlacklistManager(ConfigFactory.parseString(configString), logger = NoLogging)

    val cache = manager.blacklistCacheFor(workflowNoGrouping)
    val rootWorkflowCache = cache.get.asInstanceOf[RootWorkflowBlacklistCache]

    rootWorkflowCache.hitCache.size() shouldBe 0
    rootWorkflowCache.bucketCache.size() shouldBe 0
  }

  it should "use root workflow level caches if no workflow-option is specified for groupings in config" in {
    val configString =
      """
        |call-caching {
        |  blacklist-cache {
        |    enabled: true
        |  }
        |}
        |""".stripMargin
    val manager = new CallCachingBlacklistManager(ConfigFactory.parseString(configString), logger = NoLogging)
    val cache = manager.blacklistCacheFor(workflowYesGrouping1)
    val _ = cache.get.asInstanceOf[RootWorkflowBlacklistCache]
  }

  it should "use root workflow level caches if no workflow-option is provided in workflow options" in {
    val configString =
      """
        |call-caching {
        |  blacklist-cache {
        |    enabled: true
        |    groupings: {
        |      workflow-option: "google_project"
        |    }
        |  }
        |}
        |""".stripMargin
    val manager = new CallCachingBlacklistManager(ConfigFactory.parseString(configString), logger = NoLogging)
    val cache = manager.blacklistCacheFor(workflowNoGrouping)
    val _ = cache.get.asInstanceOf[RootWorkflowBlacklistCache]
  }

  it should "use a grouping cache if there is a workflow-option in config and its value exists as a key in workflow options" in {
    val configString =
      """
        |call-caching {
        |  blacklist-cache {
        |    enabled: true
        |    groupings: {
        |      workflow-option: "google_project"
        |    }
        |  }
        |}
        |""".stripMargin
    val manager = new CallCachingBlacklistManager(ConfigFactory.parseString(configString), logger = NoLogging)
    val cache1 = manager.blacklistCacheFor(workflowYesGrouping1).get
    val _ = cache1.asInstanceOf[GroupingBlacklistCache]

    val cache2 = manager.blacklistCacheFor(workflowYesGrouping2).get
    System.identityHashCode(cache1) shouldEqual System.identityHashCode(cache2)
  }
}
