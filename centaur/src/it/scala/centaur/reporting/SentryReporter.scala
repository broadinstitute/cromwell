package centaur.reporting

import cats.effect.IO
import centaur.reporting.SentryReporter._
import centaur.test.CentaurTestException
import com.fasterxml.jackson.core.`type`.TypeReference
import com.fasterxml.jackson.databind.ObjectMapper
import com.typesafe.config.Config
import io.sentry.dsn.Dsn
import io.sentry.event.helper.ContextBuilderHelper
import io.sentry.{DefaultSentryClientFactory, SentryClient}

import scala.util.Try

/**
  * Logs errors to a sentry dsn read from the passed in config.
  *
  * The dsn is NOT loaded using normal sentry lookup procedures.
  */
class SentryReporter(override val name: String, config: Config) extends ErrorReporter {

  val dsn: Dsn = new Dsn(config.getString("dsn"))

  override lazy val destination: String = {
    val portInfo: String = {
      Option(dsn.getPort).filterNot(_ == -1) flatMap {
        case 80 if dsn.getProtocol == "http" => None
        case 443 if dsn.getProtocol == "https" => None
        case port => Option(port)
      } map (":" + _) getOrElse ""
    }

    s"${dsn.getProtocol}://${dsn.getHost}$portInfo"
  }

  override def logCentaurFailure(testEnvironment: TestEnvironment,
                                 ciEnvironment: CiEnvironment,
                                 centaurTestException: CentaurTestException): IO[Unit] = {
    withSentryClient { sentryClient =>
      addTestEnvironment(sentryClient, testEnvironment)
      addCiEnvironment(sentryClient, ciEnvironment)
      addCentaurTestException(sentryClient, centaurTestException)
      sentryClient.sendException(centaurTestException)
    }
  }

  private def withSentryClient[A](block: SentryClient => A): IO[A] = {
    val sentryClientIo = IO(sentryReporterClientFactory.createSentryClient(dsn))
    sentryClientIo.bracket {
      sentryClient =>
        IO {
          block(sentryClient)
        }
    } {
      sentryClient =>
        IO {
          Try(sentryClient.closeConnection())
          ()
        }
    }
  }
}

object SentryReporter {
  private val CiEnvironmentPrefix = "ci_env_"

  /** Adds relevant tags for the test environment. */
  private def addTestEnvironment(sentryClient: SentryClient, testEnvironment: TestEnvironment): Unit = {
    sentryClient.addTag(CiEnvironmentPrefix + "test_name", testEnvironment.testName)
    sentryClient.addTag(CiEnvironmentPrefix + "attempt", String.valueOf(testEnvironment.attempt + 1))
  }

  /** Adds relevant tags for the CI environment, when present. */
  private def addCiEnvironment(sentryClient: SentryClient, ciEnvironment: CiEnvironment): Unit = {
    sentryClient.setEnvironment(ciEnvironment.provider.getOrElse("unknown"))
    addEnvironmentTag(sentryClient, "is_ci", ciEnvironment.isCi)
    addEnvironmentTag(sentryClient, "is_cron", ciEnvironment.isCron)
    addEnvironmentTag(sentryClient, "is_secure", ciEnvironment.isSecure)
    addEnvironmentTag(sentryClient, "type", ciEnvironment.`type`)
    addEnvironmentTag(sentryClient, "branch", ciEnvironment.branch)
    addEnvironmentTag(sentryClient, "event", ciEnvironment.event)
    addEnvironmentTag(sentryClient, "tag", ciEnvironment.tag)
    addEnvironmentTag(sentryClient, "number", ciEnvironment.number)
    addEnvironmentTag(sentryClient, "provider", ciEnvironment.provider)
    addEnvironmentTag(sentryClient, "os", ciEnvironment.os)
    addEnvironmentTag(sentryClient, "url", ciEnvironment.url)
  }

  /** Adds a tag for value, when present. */
  private def addEnvironmentTag[A](sentryClient: SentryClient, key: String, valueOption: Option[A]): Unit = {
    valueOption.map(_.toString.trim).filterNot(_.isEmpty) foreach { value =>
      sentryClient.addTag(CiEnvironmentPrefix + key, value)
    }
  }

  /** Add the workflow ID and the metadata json to the sentry report. */
  private def addCentaurTestException(sentryClient: SentryClient, centaurTestException: CentaurTestException): Unit = {
    centaurTestException.workflowIdOption foreach { workflowId =>
      sentryClient.addExtra("workflow_id", workflowId)
    }

    centaurTestException.metadataJsonOption foreach { metadataJson =>
      val mapTypeReference = new TypeReference[java.util.Map[String, AnyRef]] {}
      val mapper = new ObjectMapper()
      // Convert the metadata to a java.util.Map compatible with sentry's addExtra
      val metadataMap = mapper.readValue[java.util.Map[String, AnyRef]](metadataJson, mapTypeReference)
      sentryClient.addExtra("metadata", metadataMap)
    }
  }

  /** Overrides the default implementation to NOT return a noop client on error. */
  private val sentryReporterClientFactory = new DefaultSentryClientFactory {
    override def createSentryClient(dsn: Dsn): SentryClient = {
      val sentryClient = new SentryClient(createConnection(dsn), getContextManager(dsn))
      sentryClient.addBuilderHelper(new ContextBuilderHelper(sentryClient))
      configureSentryClient(sentryClient, dsn)
    }
  }

}
