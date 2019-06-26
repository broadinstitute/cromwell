package cromwell.services.healthmonitor.impl.workbench

import java.net.URL

import akka.actor.ActorRef
import cats.data.Validated.{Invalid, Valid}
import cats.instances.future._
import cats.syntax.functor._
import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.gax.retrying.RetrySettings
import com.google.api.services.genomics.v2alpha1.GenomicsScopes
import com.google.api.services.storage.StorageScopes
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor.{MonitoredSubsystem, OkStatus, SubsystemStatus}
import cromwell.services.healthmonitor.impl.common.{DockerHubMonitor, EngineDatabaseMonitor}
import cromwell.services.healthmonitor.impl.workbench.WorkbenchHealthMonitorServiceActor._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

/**
  * A health monitor implementation which will monitor external services used in the Workbench ecosystem, such
  * as GCS and PAPI. This implementation makes some assumptions of Cromwell's configuration which will be true
  * in a Workbench scenario but YMMV otherwise. Caveat emptor and all of that fun stuff.
  */
abstract class WorkbenchHealthMonitorServiceActor(val serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
  extends ProtoHealthMonitorServiceActor
    with DockerHubMonitor
    with EngineDatabaseMonitor {
  override implicit val system = context.system

  private lazy val papiBackendConfigurations = serviceConfig.as[Set[String]]("check-papi-backends").map(WorkbenchHealthMonitorServiceActor.PapiConfiguration.fromBackendNameValue(_, serviceConfig, globalConfig))

  def papiMonitoredSubsystem(papiConfiguration: PapiConfiguration): MonitoredSubsystem = {
    MonitoredSubsystem(papiConfiguration.backendName, () => checkPapi(papiConfiguration))
  }

  protected lazy val Gcs = MonitoredSubsystem("GCS", () => checkGcs())
  protected lazy val PapiSubsystems = papiBackendConfigurations map papiMonitoredSubsystem

  lazy val googleConfig = GoogleConfiguration(globalConfig)

  lazy val googleAuthName = serviceConfig.as[Option[String]]("google-auth-name").getOrElse("application-default")

  lazy val googleAuth = googleConfig.auth(googleAuthName) match {
    case Valid(a) => a
    case Invalid(e) => throw new IllegalArgumentException("Unable to configure WorkbenchHealthMonitor: " + e.toList.mkString(", "))
  }

  /**
    * Demonstrates connectivity to GCS by stat-ing a bucket
    */
  private def checkGcs(): Future[SubsystemStatus] = {
    // For any expected production usage of this check, the GCS bucket should be public read */
    val gcsBucketToCheck = serviceConfig.as[String]("gcs-bucket-to-check")
    val storageScopes = List(StorageScopes.DEVSTORAGE_READ_ONLY)
    val storage = Future(googleAuth.credentials(storageScopes)) map { credentials =>
      GcsStorage.gcsStorage(googleConfig.applicationName, credentials, RetrySettings.newBuilder().build())
    }
    storage map { _.buckets.get(gcsBucketToCheck).execute() } as OkStatus
  }

  private def checkPapi(papiConfiguration: PapiConfiguration): Future[SubsystemStatus] = {
    val papiConfig = papiConfiguration.papiConfig
    val papiProviderConfig = papiConfiguration.papiProviderConfig

    val endpointUrl = new URL(papiConfig.as[String]("genomics.endpoint-url"))
    val papiProjectId = papiConfig.as[String]("project")

    val check = for {
      credentials <- Future(googleAuth.credentials(List(GenomicsScopes.GENOMICS)))
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
