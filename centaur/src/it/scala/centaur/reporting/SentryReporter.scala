package centaur.reporting

import cats.effect.IO
import centaur.reporting.SentryReporter._
import centaur.test.CentaurTestException
import centaur.test.metadata.CallAttemptFailure
import com.fasterxml.jackson.core.`type`.TypeReference
import com.fasterxml.jackson.databind.ObjectMapper
import cromwell.core.WorkflowMetadataKeys
import io.sentry.dsn.Dsn
import io.sentry.event.helper.ContextBuilderHelper
import io.sentry.event.interfaces.ExceptionInterface
import io.sentry.event.{Event, EventBuilder}
import io.sentry.{DefaultSentryClientFactory, SentryClient}

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext
import scala.util.Try

/**
  * Logs errors to a sentry dsn read from the passed in config.
  *
  * The dsn is NOT loaded using normal sentry lookup procedures.
  */
class SentryReporter(override val params: ErrorReporterParams) extends ErrorReporter {

  val dsn: Dsn = new Dsn(params.reporterConfig.getString("dsn"))

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
                                 centaurTestException: CentaurTestException)
                                (implicit executionContext: ExecutionContext): IO[Unit] = {
    for {
      callAttemptFailures <- CallAttemptFailure.buildFailures(centaurTestException.metadataJsonOption)
      _ <- sendSentryFailure(testEnvironment, ciEnvironment, centaurTestException, callAttemptFailures)
    } yield ()
  }

  private def sendSentryFailure(testEnvironment: TestEnvironment,
                                ciEnvironment: CiEnvironment,
                                centaurTestException: CentaurTestException,
                                callAttemptFailures: Vector[CallAttemptFailure]): IO[Unit] = {
    withSentryClient { sentryClient =>
      addTestEnvironment(sentryClient, testEnvironment)
      addCiEnvironment(sentryClient, ciEnvironment)
      addCentaurTestException(sentryClient, centaurTestException)
      addCallAttemptFailures(sentryClient, callAttemptFailures)

      // Updated version of
      // https://github.com/getsentry/sentry-java/blob/v1.7.5/sentry/src/main/java/io/sentry/SentryClient.java#L200-L203
      val eventBuilder = new EventBuilder()
        .withMessage(centaurTestException.getMessage)
        .withLevel(Event.Level.ERROR)
        .withSentryInterface(new ExceptionInterface(centaurTestException))
        // https://docs.sentry.io/learn/rollups/#custom-grouping
        .withFingerprint("{{ default }}", providerOrDefault(ciEnvironment))

      sentryClient.sendEvent(eventBuilder)
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
  private val CiDefaultProvider = "unknown"

  private val HighlightedMetadataKeys = List(
    WorkflowMetadataKeys.SubmissionTime,
    WorkflowMetadataKeys.StartTime,
    WorkflowMetadataKeys.EndTime,
    WorkflowMetadataKeys.Status
  )

  private val HighlightedSubmissionKeys = List(
    WorkflowMetadataKeys.SubmissionSection_Workflow,
    WorkflowMetadataKeys.SubmissionSection_WorkflowUrl,
    WorkflowMetadataKeys.SubmissionSection_Root,
    WorkflowMetadataKeys.SubmissionSection_Inputs,
    WorkflowMetadataKeys.SubmissionSection_Options,
    WorkflowMetadataKeys.SubmissionSection_Imports,
    WorkflowMetadataKeys.SubmissionSection_WorkflowType,
    WorkflowMetadataKeys.SubmissionSection_Labels,
    WorkflowMetadataKeys.SubmissionSection_WorkflowTypeVersion
  )

  private def providerOrDefault(ciEnvironment: CiEnvironment): String = {
    ciEnvironment.provider.getOrElse(CiDefaultProvider)
  }

  /** Adds relevant tags for the test environment. */
  private def addTestEnvironment(sentryClient: SentryClient, testEnvironment: TestEnvironment): Unit = {
    sentryClient.addTag(CiEnvironmentPrefix + "test_name", testEnvironment.name)
    sentryClient.addTag(CiEnvironmentPrefix + "test_attempt", String.valueOf(testEnvironment.attempt + 1))
  }

  /** Adds relevant tags for the CI environment, when present. */
  private def addCiEnvironment(sentryClient: SentryClient, ciEnvironment: CiEnvironment): Unit = {
    sentryClient.setEnvironment(providerOrDefault(ciEnvironment))
    addEnvironmentTag(sentryClient, "is_ci", ciEnvironment.isCi)
    addEnvironmentTag(sentryClient, "is_cron", ciEnvironment.isCron)
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
      // Note: Metadata will be truncated in Sentry.
      // https://github.com/getsentry/sentry-java/blob/v1.7.5/sentry/src/main/java/io/sentry/marshaller/json/SentryJsonGenerator.java#L27-L30
      sentryClient.addExtra("metadata", metadataMap)

      metadataMap.asScala foreach {
        case (key, value) if HighlightedMetadataKeys.contains(key) =>
          sentryClient.addExtra(key, value)
        case (WorkflowMetadataKeys.SubmissionSection, value) =>
          value.asInstanceOf[java.util.Map[String, AnyRef]].asScala foreach {
            case (submissionKey, submissionValue) if HighlightedSubmissionKeys.contains(submissionKey) =>
              sentryClient.addExtra(WorkflowMetadataKeys.SubmissionSection + "." + submissionKey, submissionValue)
          }
        case _ => /* ignore */
      }
    }
  }

  private def addCallAttemptFailures(sentryClient: SentryClient,
                                     callAttemptFailures: Vector[CallAttemptFailure]): Unit = {
    if (callAttemptFailures.nonEmpty) {
      val failures = callAttemptFailures map { failure =>
        Map(
          "call_fully_qualified_name" -> Option(failure.callFullyQualifiedName),
          "job_index" -> Option(failure.jobIndex),
          "job_attempt" -> Option(failure.jobAttempt),
          "message" -> Option(failure.message),
          "start" -> failure.startOption,
          "end" -> failure.endOption,
          "stdout" -> failure.stdoutOption,
          "stderr" -> failure.stderrOption,
          "call_root" -> failure.callRootOption,
        ) collect {
          case (key, Some(value)) => (key, value)
        }
      }
      sentryClient.addExtra(s"call_failures", failures.map(_.asJava).toArray)
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
