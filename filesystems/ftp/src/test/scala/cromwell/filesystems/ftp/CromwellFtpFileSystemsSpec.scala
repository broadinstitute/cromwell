package cromwell.filesystems.ftp

import cloud.nio.impl.ftp.FtpFileSystemsConfiguration.Active
import com.typesafe.config.ConfigFactory
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._

class CromwellFtpFileSystemsSpec extends AnyFlatSpec with Matchers {

  behavior of "CromwellFtpFileSystemsSpec"

  it should "parse configuration" in {
    val config = ConfigFactory.parseString(
      """cache-ttl = 10 days
        |obtain-connection-timeout = 12 hours
        |max-connection-per-server-per-user = 1
        |idle-connection-timeout = 14 hours
        |connection-port: 212
        |connection-mode = "active" """.stripMargin)
    
    val fs = new CromwellFtpFileSystems(config)
    fs.ftpFileSystems.config.cacheTTL shouldBe 10.days
    fs.ftpFileSystems.config.leaseTimeout shouldBe Some(12.hours)
    // 1 should be upped to 2
    fs.ftpFileSystems.config.capacity shouldBe 2
    fs.ftpFileSystems.config.idleConnectionTimeout shouldBe 14.hours
    fs.ftpFileSystems.config.connectionPort shouldBe 212
    fs.ftpFileSystems.config.connectionMode shouldBe Active
  }

}
