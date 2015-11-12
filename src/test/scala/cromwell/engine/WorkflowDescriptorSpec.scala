package cromwell.engine

import com.typesafe.config.ConfigFactory
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._

class WorkflowDescriptorSpec extends FlatSpec with Matchers {

  val dummyWdl = "workflow w {}"
  val configMap = Map("backend.backend" -> "local")

  it should "honor configuration and workflow options for call-caching" in {
    val callCachingOnTable = Table(
      ("config", "wfOptions"),

      (configMap, "{\"use-cache\": true}"),
      (configMap + ("call-caching.enabled" -> "true"), "{}"),
      (configMap + ("call-caching.enabled" -> "true"), "{\"use-cache\": true}"),
      (configMap + ("call-caching.enabled" -> "false"), "{\"use-cache\": true}")
    )

    val callCachingOffTable = Table(
      ("config", "wfOptions"),

      (configMap, "{}"),
      (configMap, "{\"use-cache\": false}"),
      (configMap + ("call-caching.enabled" -> "false"), "{}"),
      (configMap + ("call-caching.enabled" -> "false"), "{\"use-cache\": false}"),
      (configMap + ("call-caching.enabled" -> "true"), "{\"use-cache\": false}")
    )

    forAll(callCachingOnTable) { (config, options) =>
      val sources = WorkflowSourceFiles(dummyWdl, "{}", options)
      val wd = new WorkflowDescriptor(WorkflowId.randomId(), sources) {
        override lazy val conf = ConfigFactory.parseMap(config.asJava)
      }
      wd.cacheCalls shouldBe true
    }

    forAll(callCachingOffTable) { (config, options) =>
      val sources = WorkflowSourceFiles(dummyWdl, "{}", options)
      val wd = new WorkflowDescriptor(WorkflowId.randomId(), sources) {
        override lazy val conf = ConfigFactory.parseMap(config.asJava)
      }
      wd.cacheCalls shouldBe false
    }
  }

}
