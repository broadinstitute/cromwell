package cromwell.filesystems.ftp

import akka.actor.ActorSystem
import cloud.nio.impl.ftp.{FtpAuthenticatedCredentials, FtpCloudNioFileSystemProvider}
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.ftp.FtpPathBuilderFactory._

import scala.concurrent.{ExecutionContext, Future}

object FtpPathBuilderFactory {
  case class FileSystemKey(host: String, configuration: FtpConfiguration)
  def credentialsFromWorkflowOptions(workflowOptions: WorkflowOptions) = {
    def getValue(key: String) = workflowOptions.get("ftp-username").toOption

    (getValue("ftp-username"), getValue("ftp-password"), getValue("ftp-account")) match {
      case (Some(username), Some(password), account) => Option(FtpAuthenticatedCredentials(username, password, account))
      case _ => None
    }
  }
}

class FtpPathBuilderFactory(globalConfig: Config, instanceConfig: Config) extends PathBuilderFactory {
  private [ftp] lazy val configFtpConfiguration = FtpConfiguration(instanceConfig)
  private lazy val defaultFtpProvider = new FtpCloudNioFileSystemProvider(instanceConfig, configFtpConfiguration.ftpCredentials)

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) = Future.successful {
    val provider = credentialsFromWorkflowOptions(options) match {
      case Some(overriddenCredentials) =>
        new FtpCloudNioFileSystemProvider(instanceConfig, overriddenCredentials)
      case _ => defaultFtpProvider
    }
    
    FtpPathBuilder(provider)
  }
}
