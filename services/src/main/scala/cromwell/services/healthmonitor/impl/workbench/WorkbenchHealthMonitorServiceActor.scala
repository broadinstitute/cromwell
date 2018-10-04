package cromwell.services.healthmonitor.impl.workbench

import java.net.URL

import akka.actor.ActorRef
import cats.data.Validated.{Invalid, Valid}
import cats.instances.future._
import cats.syntax.functor._
import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.gax.retrying.RetrySettings
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.services.healthmonitor.HealthMonitorServiceActor
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{MonitoredSubsystem, OkStatus, SubsystemStatus}
import cromwell.services.healthmonitor.impl.common.{DockerHubMonitor, EngineDatabaseMonitor}
import cromwell.services.healthmonitor.impl.workbench.WorkbenchHealthMonitorServiceActor._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

/**
  * A health monitor implementation which will monitor external services used in the Workbench ecosystem, such
  * as GCS and PAPI. This implementation makes some assumptions of Cromwell's configuration which will be true
  * in a Workbench scenario but YMMV otherwise. Caveat emptor and all of that fun stuff.
  */
class WorkbenchHealthMonitorServiceActor(val serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
  extends HealthMonitorServiceActor
    with DockerHubMonitor
    with EngineDatabaseMonitor {
  override implicit val system = context.system

  private lazy val papiV1ConfigOption =
    PapiConfiguration.fromBackendNameKey("papi-v1-backend-name", serviceConfig, globalConfig)
  private lazy val papiV2ConfigOption =
    PapiConfiguration.fromBackendNameKey("papi-v2-backend-name", serviceConfig, globalConfig)
  private lazy val papiConfigOption =
    PapiConfiguration.fromBackendNameKey("papi-backend-name", serviceConfig, globalConfig)
  private lazy val defaultJesConfigOption = (papiConfigOption, papiV1ConfigOption, papiV2ConfigOption) match {
    case (None, None, None) => Option(PapiConfiguration.fromBackendNameValue("Jes", serviceConfig, globalConfig))
    case _ => None
  }

  def papiMonitoredSubsystem(name: String)(papiConfiguration: PapiConfiguration): MonitoredSubsystem = {
    MonitoredSubsystem(name, () => checkPapi(papiConfiguration))
  }

  private lazy val Gcs = MonitoredSubsystem("GCS", () => checkGcs())
  private lazy val PapiV1Option = papiV1ConfigOption map papiMonitoredSubsystem("PAPIV1")
  private lazy val PapiV2Option = papiV2ConfigOption map papiMonitoredSubsystem("PAPIV2")
  private lazy val PapiOption = papiConfigOption map papiMonitoredSubsystem("PAPI")
  private lazy val DefaultJesOption = defaultJesConfigOption map papiMonitoredSubsystem("PAPI")

  override lazy val subsystems: Set[MonitoredSubsystem] = Set(
    Option(DockerHub),
    Option(EngineDb),
    PapiOption,
    PapiV1Option,
    PapiV2Option,
    DefaultJesOption,
    Option(Gcs)
  ) collect {
    case Some(value) => value
  }

  val googleConfig = GoogleConfiguration(globalConfig)

  val googleAuthName = serviceConfig.as[Option[String]]("google-auth-name").getOrElse("application-default")

  val googleAuth = googleConfig.auth(googleAuthName) match {
    case Valid(a) => a
    case Invalid(e) => throw new IllegalArgumentException("Unable to configure WorkbenchHealthMonitor: " + e.toList.mkString(", "))
  }

  /**
    * Demonstrates connectivity to GCS by stat-ing a bucket
    */
  private def checkGcs(): Future[SubsystemStatus] = {
    // For any expected production usage of this check, the GCS bucket should be public read */
    val gcsBucketToCheck = serviceConfig.as[String]("gcs-bucket-to-check")
    val storage = Future(googleAuth.pipelinesApiCredentials(GoogleAuthMode.NoOptionLookup)) map { credentials =>
      GcsStorage.gcsStorage(googleConfig.applicationName, credentials, RetrySettings.newBuilder().build())
    }
    storage map { _.buckets.get(gcsBucketToCheck).execute() } as OkStatus
  }

  private def checkPapi(papiConfiguration: PapiConfiguration): Future[SubsystemStatus] = {
    checkPapi(papiConfiguration.papiConfig, papiConfiguration.papiProviderConfig)
  }

  /**
    * Demonstrates connectivity to Google Pipelines API (PAPI) by making sure it can access an authenticated endpoint
    */
  private def checkPapi(papiConfig: Config, papiProviderConfig: Config): Future[SubsystemStatus] = {
    val endpointUrl = new URL(papiConfig.as[String]("genomics.endpoint-url"))
    val papiProjectId = papiConfig.as[String]("project")

    val check = for {
      credentials <- Future(googleAuth.pipelinesApiCredentials(GoogleAuthMode.NoOptionLookup))
      genomicsChecker = if (papiProviderConfig.as[String]("actor-factory").contains("v2alpha1"))
        GenomicsCheckerV2(googleConfig.applicationName, googleAuth, endpointUrl, credentials, papiProjectId)
      else
        GenomicsCheckerV1(googleConfig.applicationName, googleAuth, endpointUrl, credentials, papiProjectId)
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
        def initialize(httpRequest: HttpRequest) = {
          delegate.initialize(httpRequest)
        }
      }
    }

    def check: Future[Unit]
  }

  case class GenomicsCheckerV1(applicationName: String,
                               authMode: GoogleAuthMode,
                               endpointUrl: URL,
                               credentials: Credentials,
                               papiProjectId: String)(implicit val ec: ExecutionContext) extends GenomicsChecker {
    val genomics = new com.google.api.services.genomics.Genomics.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      httpInitializer(credentials))
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def check = Future {
      genomics.pipelines().list().setProjectId(papiProjectId).setPageSize(1).execute()
      ()
    }
  }

  case class GenomicsCheckerV2(applicationName: String,
                               authMode: GoogleAuthMode,
                               endpointUrl: URL,
                               credentials: Credentials,
                               papiProjectId: String)(implicit val ec: ExecutionContext) extends GenomicsChecker {
    val genomics = new com.google.api.services.genomics.v2alpha1.Genomics.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      httpInitializer(credentials))
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build

    override def check = Future {
      // https://cloud.google.com/genomics/reference/rest/#rest-resource-v2alpha1projectsoperations
      genomics.projects().operations().list(s"projects/$papiProjectId/operations").setPageSize(1).execute()
      ()
    }
  }

  case class PapiConfiguration(backendName: String, papiConfig: Config, papiProviderConfig: Config)

  object PapiConfiguration {
    def fromBackendNameKey(backendNameKey: String,
                           serviceConfig: Config,
                           globalConfig: Config): Option[PapiConfiguration] = {
      serviceConfig.as[Option[String]](backendNameKey) map {
        fromBackendNameValue(_, serviceConfig, globalConfig)
      }
    }

    def fromBackendNameValue(papiBackendName: String,
                             serviceConfig: Config,
                             globalConfig: Config): PapiConfiguration = {
      val papiProviderConfig: Config = globalConfig.as[Config](s"backend.providers.$papiBackendName")
      val papiConfig: Config = papiProviderConfig.as[Config]("config")
      PapiConfiguration(papiBackendName, papiConfig, papiProviderConfig)
    }
  }
}
