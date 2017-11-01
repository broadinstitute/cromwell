package cromwell.core

import java.net.URL

import com.typesafe.config.{Config, ConfigFactory}
import org.scalatest.{FlatSpecLike, Matchers, Suite}

import scala.collection.JavaConverters._

class ReferenceConfSpec extends Suite with FlatSpecLike with Matchers {

  behavior of "reference.conf"

  it should "check for colliding cromwell keys" in {
    val classLoader = Thread.currentThread().getContextClassLoader
    val references = classLoader.getResources("reference.conf").asScala.toSeq
    val (fileUrls, nonFileUrls) = references.partition(_.getProtocol == "file")

    if (fileUrls.isEmpty)
      fail(s"Could not locate any file urls within:\n  ${nonFileUrls.mkString("\n  ")}")

    for (fileUrl <- fileUrls) {
      val fileKeys = getKeys(fileUrl)
      for (nonFileUrl <- nonFileUrls) {
        val nonFileKeys = getKeys(nonFileUrl)

        /*
        Key definitions go in reference.conf.
        Key overrides go in application.conf.
        FYI: Stanzas, for ex: "my { empty { stanza {} } }" never collide.
        FYI: Lists count as values, not stanzas, and will collide.

        We should never define keys within cromwell's reference.conf that collide with keys from a library
        reference.conf. Overriding a key's value should be placed into application.conf.

        Examples:
        reference.conf definitions:
          akka.dispatchers.my-new-dispatcher.type = Dispatcher
          akka.dispatchers.my-new-dispatcher.executor = "fork-join-executor"

        application.conf overrides:
          akka.loggers = ["my.list", "of.loggers"]
        */
        val overlappingKeys = fileKeys.intersect(nonFileKeys).toSeq.sorted
        if (overlappingKeys.nonEmpty)
          fail(
            s"""|Key(s) overlapping
                |=> ${overlappingKeys.mkString("=> ")}
                |between $fileUrl
                |    and $nonFileUrl
                |""".stripMargin)
      }
    }
  }

  private def getKeys(configUrl: URL): Set[String] = {
    getKeys(ConfigFactory.parseURL(configUrl))
  }

  private def getKeys(config: Config): Set[String] = {
    config.entrySet().asScala.map(_.getKey).toSet
  }
}
