package cloud.nio.impl.drs

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

final case class DrsConfig(marthaUrl: String,
                           numRetries: Int,
                           waitInitial: FiniteDuration,
                           waitMaximum: FiniteDuration,
                           waitMultiplier: Double,
                           waitRandomizationFactor: Double,
                          )

object DrsConfig {
  // If you update these values also update Filesystems.md!
  private val DefaultNumRetries = 3
  private val DefaultWaitInitial = 30.seconds
  private val DefaultWaitMaximum = 5.minutes
  private val DefaultWaitMultiplier = 2.0d
  private val DefaultWaitRandomizationFactor = 0.1

  private val EnvMarthaUrl = "MARTHA_URL"
  private val EnvMarthaNumRetries = "MARTHA_NUM_RETRIES"
  private val EnvMarthaWaitInitialSeconds = "MARTHA_WAIT_INITIAL_SECONDS"
  private val EnvMarthaWaitMaximumSeconds = "MARTHA_WAIT_MAXIMUM_SECONDS"
  private val EnvMarthaWaitMultiplier = "MARTHA_WAIT_MULTIPLIER"
  private val EnvMarthaWaitRandomizationFactor = "MARTHA_WAIT_RANDOMIZATION_FACTOR"

  def fromConfig(marthaConfig: Config): DrsConfig = {
    DrsConfig(
      marthaUrl = marthaConfig.getString("url"),
      numRetries = marthaConfig.getOrElse("num-retries", DefaultNumRetries),
      waitInitial = marthaConfig.getOrElse("wait-initial", DefaultWaitInitial),
      waitMaximum = marthaConfig.getOrElse("wait-maximum", DefaultWaitMaximum),
      waitMultiplier = marthaConfig.getOrElse("wait-multiplier", DefaultWaitMultiplier),
      waitRandomizationFactor =
        marthaConfig.getOrElse("wait-randomization-factor", DefaultWaitRandomizationFactor),
    )
  }

  def fromEnv(env: Map[String, String]): DrsConfig = {
    DrsConfig(
      marthaUrl = env(EnvMarthaUrl),
      numRetries = env.get(EnvMarthaNumRetries).map(_.toInt).getOrElse(DefaultNumRetries),
      waitInitial = env.get(EnvMarthaWaitInitialSeconds).map(_.toLong.seconds).getOrElse(DefaultWaitInitial),
      waitMaximum = env.get(EnvMarthaWaitMaximumSeconds).map(_.toLong.seconds).getOrElse(DefaultWaitMaximum),
      waitMultiplier = env.get(EnvMarthaWaitMultiplier).map(_.toDouble).getOrElse(DefaultWaitMultiplier),
      waitRandomizationFactor =
        env.get(EnvMarthaWaitRandomizationFactor).map(_.toDouble).getOrElse(DefaultWaitRandomizationFactor),
    )
  }

  def toEnv(drsConfig: DrsConfig): Map[String, String] = {
    Map(
      EnvMarthaUrl -> drsConfig.marthaUrl,
      EnvMarthaNumRetries -> s"${drsConfig.numRetries}",
      EnvMarthaWaitInitialSeconds -> s"${drsConfig.waitInitial.toSeconds}",
      EnvMarthaWaitMaximumSeconds -> s"${drsConfig.waitMaximum.toSeconds}",
      EnvMarthaWaitMultiplier -> s"${drsConfig.waitMultiplier}",
      EnvMarthaWaitRandomizationFactor -> s"${drsConfig.waitRandomizationFactor}",
    )
  }
}
