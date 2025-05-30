package cromwell.services.healthmonitor.impl.workbench

import java.net.URL

import akka.actor.ActorRef
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.functor._
import cats.instances.future._
import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.lifesciences.v2beta.CloudLifeSciencesScopes
import com.google.api.services.lifesciences.v2beta.CloudLifeSciences
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor.{MonitoredSubsystem, OkStatus, SubsystemStatus}
import cromwell.services.healthmonitor.impl.common.EngineDatabaseMonitor
import cromwell.services.healthmonitor.impl.workbench.WorkbenchHealthMonitorServiceActor._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

/**
  * A health monitor implementation which will monitor external services used in the Workbench ecosystem, such
  * as GCS and PAPI. This implementation makes some assumptions of Cromwell's configuration which will be true
  * in a Workbench scenario but YMMV otherwise. Caveat emptor and all of that fun stuff.
  */
abstract class WorkbenchHealthMonitorServiceActor(val serviceConfig: Config,
                                                  globalConfig: Config,
                                                  serviceRegistryActor: ActorRef
) extends ProtoHealthMonitorServiceActor
    with EngineDatabaseMonitor {
  private lazy val papiBackendConfigurations = serviceConfig
    .as[Set[String]]("check-papi-backends")
    .map(WorkbenchHealthMonitorServiceActor.PapiConfiguration.fromBackendNameValue(_, serviceConfig, globalConfig))

  def papiMonitoredSubsystem(papiConfiguration: PapiConfiguration): MonitoredSubsystem =
    MonitoredSubsystem(papiConfiguration.backendName, () => checkPapi(papiConfiguration))

  protected lazy val PapiSubsystems = papiBackendConfigurations map papiMonitoredSubsystem

  lazy val googleConfig = GoogleConfiguration(globalConfig)

  lazy val googleAuthName = serviceConfig.as[Option[String]]("google-auth-name").getOrElse("application-default")
  lazy val googleAuth = getGoogleAuthConfigurationOrFail(googleAuthName)

  private def getGoogleAuthConfigurationOrFail(googleAuthName: String): GoogleAuthMode =
    googleConfig.auth(googleAuthName) match {
      case Valid(a) => a
      case Invalid(e) =>
        throw new IllegalArgumentException("Unable to configure WorkbenchHealthMonitor: " + e.toList.mkString(", "))
    }

  private def checkPapi(papiConfiguration: PapiConfiguration): Future[SubsystemStatus] = {
    val papiConfig = papiConfiguration.papiConfig
    val endpointUrl = new URL(papiConfig.as[String]("genomics.endpoint-url"))
    val papiProjectId = papiConfig.as[String]("project")
    val location = papiConfig.as[String]("genomics.location")
    val check = for {
      credentials <- Future(googleAuth.credentials(List(CloudLifeSciencesScopes.CLOUD_PLATFORM)))
      genomicsChecker = GenomicsCheckerV2Beta(googleConfig.applicationName,
                                              googleAuth,
                                              endpointUrl,
                                              location,
                                              credentials,
                                              papiProjectId
      )
      checked <- genomicsChecker.check
    } yield checked

    check as OkStatus
  }
}

object WorkbenchHealthMonitorServiceActor {
  sealed trait GenomicsChecker {
    protected def httpInitializer(credentials: Credentials) = {
      val delegate = new HttpCredentialsAdapter(credentials)
      new HttpRequestInitializer() {
        def initialize(httpRequest: HttpRequest) =
          delegate.initialize(httpRequest)
      }
    }

    def check: Future[Unit]
  }

  case class GenomicsCheckerV2Beta(applicationName: String,
                                   authMode: GoogleAuthMode,
                                   endpointUrl: URL,
                                   location: String,
                                   credentials: Credentials,
                                   papiProjectId: String
  )(implicit val ec: ExecutionContext)
      extends GenomicsChecker {
    val lifeSciences = new CloudLifeSciences.Builder(GoogleAuthMode.httpTransport,
                                                     GoogleAuthMode.jsonFactory,
                                                     httpInitializer(credentials)
    )
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def check = Future {
      // https://cloud.google.com/life-sciences/docs/reference/rest/v2beta/projects.locations.operations
      lifeSciences
        .projects()
        .locations()
        .operations()
        .list(s"projects/$papiProjectId/locations/$location")
        .setPageSize(1)
        .execute()
      ()
    }
  }

  case class PapiConfiguration(backendName: String, papiConfig: Config, papiProviderConfig: Config)

  object PapiConfiguration {
    def fromBackendNameKey(backendNameKey: String,
                           serviceConfig: Config,
                           globalConfig: Config
    ): Option[PapiConfiguration] =
      serviceConfig.as[Option[String]](backendNameKey) map {
        fromBackendNameValue(_, serviceConfig, globalConfig)
      }

    def fromBackendNameValue(papiBackendName: String,
                             serviceConfig: Config,
                             globalConfig: Config
    ): PapiConfiguration = {
      val papiProviderConfig: Config = globalConfig.as[Config](s"backend.providers.$papiBackendName")
      val papiConfig: Config = papiProviderConfig.as[Config]("config")
      PapiConfiguration(papiBackendName, papiConfig, papiProviderConfig)
    }
  }
}
