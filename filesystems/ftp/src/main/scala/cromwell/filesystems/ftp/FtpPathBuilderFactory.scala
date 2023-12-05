package cromwell.filesystems.ftp

import akka.actor.ActorSystem
import cloud.nio.impl.ftp.{FtpAuthenticatedCredentials, FtpCloudNioFileSystemProvider}
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import cromwell.filesystems.ftp.FtpPathBuilderFactory._

import scala.concurrent.{ExecutionContext, Future}

object FtpPathBuilderFactory {
  object WorkflowOptions {
    val FtpUsername = "ftp-username"
    val FtpPassword = "ftp-password"
    val FtpAccount = "ftp-account"
  }

  def credentialsFromWorkflowOptions(workflowOptions: WorkflowOptions) = {
    def getValue(key: String) = workflowOptions.get(key).toOption

    (getValue(WorkflowOptions.FtpUsername),
     getValue(WorkflowOptions.FtpPassword),
     getValue(WorkflowOptions.FtpAccount)
    ) match {
      case (Some(username), Some(password), account) => Option(FtpAuthenticatedCredentials(username, password, account))
      case _ => None
    }
  }
}

class FtpPathBuilderFactory(globalConfig: Config,
                            instanceConfig: Config,
                            cromwellFtpFileSystems: CromwellFtpFileSystems
) extends PathBuilderFactory {
  private[ftp] lazy val configFtpConfiguration = FtpInstanceConfiguration(instanceConfig)
  private lazy val defaultFtpProvider = new FtpCloudNioFileSystemProvider(instanceConfig,
                                                                          configFtpConfiguration.ftpCredentials,
                                                                          cromwellFtpFileSystems.ftpFileSystems
  )

  override def withOptions(options: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext) =
    Future.successful {
      val provider = credentialsFromWorkflowOptions(options) match {
        case Some(overriddenCredentials) =>
          new FtpCloudNioFileSystemProvider(instanceConfig,
                                            overriddenCredentials,
                                            cromwellFtpFileSystems.ftpFileSystems
          )
        case _ => defaultFtpProvider
      }

      FtpPathBuilder(provider)
    }
}
