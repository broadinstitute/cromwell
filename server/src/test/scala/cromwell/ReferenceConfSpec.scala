package cromwell

import java.net.URL

import com.typesafe.config.ConfigFactory
import org.scalatest.{FlatSpecLike, Matchers, Suite}

import scala.collection.JavaConverters._

class ReferenceConfSpec extends Suite with FlatSpecLike with Matchers {

  behavior of "reference.conf"

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
  it should "check for colliding cromwell keys" in {
    val (fileUrls, nonFileUrls) = partitionConfigFiles("reference.conf")

    for (fileUrl <- fileUrls) {
      val fileKeys = getKeys(fileUrl)
      for (nonFileUrl <- nonFileUrls) {
        val nonFileKeys = getKeys(nonFileUrl)

        val overlappingKeys = fileKeys.intersect(nonFileKeys).toSeq.sorted
        if (overlappingKeys.nonEmpty) {
          fail(
            s"""|Key(s) overlapping
                |=> ${overlappingKeys.mkString("=> ")}
                |between $fileUrl
                |    and $nonFileUrl
                |""".stripMargin)
        }
      }
    }
  }

  /*
  New keys should be defined in reference.conf.
  Overrides of existing keys should go in application.conf.

  We should never define new root keys within cromwell's application.conf.

  Bad:
    reference.conf does not do definition:
      # Does not contain a service stanza

    application.conf defining a new stanza:
      services {
        # defines new services...
      }

  Good:
    reference.conf definitions:
      services {
        # defines new services...
      }

    application.conf other overrides:
      # No new cromwell definitions, but does override reference values from libraries...
   */
  it should "check for new root keys in application.conf" in {
    val applicationFileRootKeys = getRootFileKeys("application.conf")
    val referenceFileRootKeys = getRootFileKeys("reference.conf")

    val newKeys = applicationFileRootKeys.diff(referenceFileRootKeys).toSeq.sorted
    if (newKeys.nonEmpty) {
      fail(s"New root key(s) only defined in application.conf: ${newKeys.mkString(", ")}")
    }
  }

  private def partitionConfigFiles(configFileName: String): (Seq[URL], Seq[URL]) = {
    val classLoader = Thread.currentThread().getContextClassLoader
    val urls = classLoader.getResources(configFileName).asScala.toSeq
    val (fileUrls, otherUrls) = urls.partition(_.getProtocol == "file")

    val mainFileUrls = fileUrls.filterNot(url => url.toExternalForm.contains("/test-classes/"))
    if (mainFileUrls.isEmpty)
      fail(s"Could not locate any main file urls within:\n  ${urls.mkString("\n  ")}")

    (mainFileUrls, otherUrls)
  }

  private def getRootFileKeys(configFileName: String): Set[String] = {
    val (fileUrls, _) = partitionConfigFiles(configFileName)
    for {
      fileUrl <- fileUrls.toSet[URL]
      rootKey <- getRootKeys(fileUrl)
    } yield rootKey
  }

  private def getRootKeys(configUrl: URL): Set[String] = {
    val config = ConfigFactory.parseURL(configUrl)
    val entries = config.root().entrySet()
    val keys = entries.asScala.map(_.getKey).toSet
    keys
  }

  private def getKeys(configUrl: URL): Set[String] = {
    val config = ConfigFactory.parseURL(configUrl)
    val entries = config.entrySet()
    val keys = entries.asScala.map(_.getKey).toSet
    keys
  }
}
