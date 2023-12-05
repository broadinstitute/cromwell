package cloud.nio.impl.drs

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps

final case class DrsConfig(drsResolverUrl: String,
                           numRetries: Int,
                           waitInitial: FiniteDuration,
                           waitMaximum: FiniteDuration,
                           waitMultiplier: Double,
                           waitRandomizationFactor: Double
)

object DrsConfig {
  // If you update these values also update Filesystems.md!
  private val DefaultNumRetries = 3
  private val DefaultWaitInitial = 30 seconds
  private val DefaultWaitMaximum = 60 seconds
  private val DefaultWaitMultiplier = 1.25d
  private val DefaultWaitRandomizationFactor = 0.1

  private val EnvDrsResolverUrl = "DRS_RESOLVER_URL"
  private val EnvDrsResolverNumRetries = "DRS_RESOLVER_NUM_RETRIES"
  private val EnvDrsResolverWaitInitialSeconds = "DRS_RESOLVER_WAIT_INITIAL_SECONDS"
  private val EnvDrsResolverWaitMaximumSeconds = "DRS_RESOLVER_WAIT_MAXIMUM_SECONDS"
  private val EnvDrsResolverWaitMultiplier = "DRS_RESOLVER_WAIT_MULTIPLIER"
  private val EnvDrsResolverWaitRandomizationFactor = "DRS_RESOLVER_WAIT_RANDOMIZATION_FACTOR"

  def fromConfig(drsResolverConfig: Config): DrsConfig =
    DrsConfig(
      drsResolverUrl = drsResolverConfig.getString("url"),
      numRetries = drsResolverConfig.getOrElse("num-retries", DefaultNumRetries),
      waitInitial = drsResolverConfig.getOrElse("wait-initial", DefaultWaitInitial),
      waitMaximum = drsResolverConfig.getOrElse("wait-maximum", DefaultWaitMaximum),
      waitMultiplier = drsResolverConfig.getOrElse("wait-multiplier", DefaultWaitMultiplier),
      waitRandomizationFactor = drsResolverConfig.getOrElse("wait-randomization-factor", DefaultWaitRandomizationFactor)
    )

  def fromEnv(env: Map[String, String]): DrsConfig =
    DrsConfig(
      drsResolverUrl = env(EnvDrsResolverUrl),
      numRetries = env.get(EnvDrsResolverNumRetries).map(_.toInt).getOrElse(DefaultNumRetries),
      waitInitial = env.get(EnvDrsResolverWaitInitialSeconds).map(_.toLong.seconds).getOrElse(DefaultWaitInitial),
      waitMaximum = env.get(EnvDrsResolverWaitMaximumSeconds).map(_.toLong.seconds).getOrElse(DefaultWaitMaximum),
      waitMultiplier = env.get(EnvDrsResolverWaitMultiplier).map(_.toDouble).getOrElse(DefaultWaitMultiplier),
      waitRandomizationFactor =
        env.get(EnvDrsResolverWaitRandomizationFactor).map(_.toDouble).getOrElse(DefaultWaitRandomizationFactor)
    )

  def toEnv(drsConfig: DrsConfig): Map[String, String] =
    Map(
      EnvDrsResolverUrl -> drsConfig.drsResolverUrl,
      EnvDrsResolverNumRetries -> s"${drsConfig.numRetries}",
      EnvDrsResolverWaitInitialSeconds -> s"${drsConfig.waitInitial.toSeconds}",
      EnvDrsResolverWaitMaximumSeconds -> s"${drsConfig.waitMaximum.toSeconds}",
      EnvDrsResolverWaitMultiplier -> s"${drsConfig.waitMultiplier}",
      EnvDrsResolverWaitRandomizationFactor -> s"${drsConfig.waitRandomizationFactor}"
    )
}
