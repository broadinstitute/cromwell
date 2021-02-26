package centaur.reporting

import cats.effect.IO
import centaur.{CentaurConfig, CromwellDatabase}
import com.typesafe.config.{Config, ConfigFactory}
import com.typesafe.scalalogging.StrictLogging
import net.ceedubs.ficus.Ficus._

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext

/**
  * A collection of reporters.
  *
  * @param rootConfig The config used to load the reporters.
  */
class ErrorReporters(rootConfig: Config) {

  private val errorReporterConfig: Config = rootConfig.getOrElse("error-reporter", ConfigFactory.empty)

  // Parameters to CromwellDatabase() should be lazy and are only initialized when needed.
  private val cromwellDatabase = new ErrorReporterCromwellDatabase(CromwellDatabase.fromConfig(CentaurConfig.conf))

  private val errorReporterNames: List[String] = {
    val providersConfig = errorReporterConfig.getOrElse("providers", ConfigFactory.empty)
    providersConfig.entrySet.asScala.map(_.getKey.split("\\.").toList.head).toList
  }

  private val errorReportersIo: IO[List[ErrorReporter]] = {
    AggregatedIo.aggregateExceptions("Errors while creating ErrorReporters", errorReporterNames.map(getErrorReporter))
  }

  val errorReporters: List[ErrorReporter] = errorReportersIo.unsafeRunSync

  /** The number of times any test should be retried. */
  val retryAttempts: Int = errorReporterConfig.getOrElse("retry-attempts", 0)

  /**
    * Log a centaur failure to this collection of reporters.
    *
    * @param testEnvironment      The test information, including the name of the test and the attemp.
    * @param ciEnvironment        The information about the CI environment running the test.
    * @param throwable            The exception that occurred while running the test.
    * @return An IO effect that will log the failure.
    */
  def logFailure(testEnvironment: TestEnvironment,
                 ciEnvironment: CiEnvironment,
                 throwable: Throwable)
                (implicit executionContext: ExecutionContext): IO[Unit] = {
    if (errorReporters.isEmpty) {
      // If the there are no reporters, then just "throw" the exception. Do not retry to run the test.
      IO.raiseError(throwable)
    } else {
      val listIo = errorReporters.map(_.logFailure(testEnvironment, ciEnvironment, throwable))
      AggregatedIo.aggregateExceptions("Errors while reporting a failure", listIo).handleErrorWith(err => {
        err.addSuppressed(throwable)
        IO.raiseError(err)
      }).void
    }
  }

  /**
    * Constructs the IO reporter by name.
    */
  private def getErrorReporter(errorReporterName: String): IO[ErrorReporter] = {
    IO {
      val clazz = errorReporterConfig.getString(s"providers.$errorReporterName.class")
      val reporterConfig = errorReporterConfig.getOrElse(s"providers.$errorReporterName.config", ConfigFactory.empty)
      val constructor = Class.forName(clazz).getConstructor(classOf[ErrorReporterParams])
      val params = ErrorReporterParams(errorReporterName, rootConfig, reporterConfig, cromwellDatabase)
      constructor.newInstance(params).asInstanceOf[ErrorReporter]
    }
  }
}

object ErrorReporters extends StrictLogging {
  private val ciEnvironment = CiEnvironment()
  private[reporting] val errorReporters = new ErrorReporters(CentaurConfig.conf)
  val retryAttempts = errorReporters.retryAttempts

  errorReporters.errorReporters foreach { errorReporter =>
    logger.info("Error reporter loaded: {} to {}", errorReporter.params.name, errorReporter.destination)
  }

  if (retryAttempts > 0)
    logger.info("Error retry count: {}", retryAttempts)

  def logFailure(testEnvironment: TestEnvironment,
                 throwable: Throwable)
                (implicit executionContext: ExecutionContext): IO[Unit] = {
    errorReporters.logFailure(testEnvironment, ciEnvironment, throwable)
  }
}
