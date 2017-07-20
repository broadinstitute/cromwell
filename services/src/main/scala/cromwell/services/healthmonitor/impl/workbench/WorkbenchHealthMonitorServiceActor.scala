package cromwell.services.healthmonitor.impl.workbench

import java.net.URL

import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.Config
import cromwell.cloudSupport.gcp.GoogleConfiguration
import cromwell.cloudSupport.gcp.gcs.GcsStorage
import cromwell.cloudSupport.gcp.genomics.GenomicsFactory
import cromwell.services.healthmonitor.HealthMonitorServiceActor
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{MonitoredSubsystem, OkStatus, SubsystemStatus}
import cromwell.services.healthmonitor.impl.CommonMonitoredSubsystems
import net.ceedubs.ficus.Ficus._

import scala.concurrent.Future

/**
  * A health monitor implementation which will monitor external services used in the Workbench ecosystem, such
  * as GCS and PAPI
  */
class WorkbenchHealthMonitorServiceActor(val serviceConfig: Config, globalConfig: Config)
  extends HealthMonitorServiceActor
    with CommonMonitoredSubsystems {
  implicit val as = context.system

  override lazy val subsystems: List[MonitoredSubsystem] = List(DockerHub, EngineDb, Papi, Gcs)

  private lazy val Gcs = MonitoredSubsystem("GCS", checkGcs _)
  private lazy val Papi = MonitoredSubsystem("PAPI", checkPapi _)

  val googleConfig = GoogleConfiguration(globalConfig)
  val papiBackendName = serviceConfig.as[Option[String]]("services.HealthMonitor.PapiBackendName").getOrElse("Jes")
  val papiConfig: Config = globalConfig.as[Config]("backend.providers.Jes.config")

  /**
    * Demonstrates connectivity to GCS by stat-ing the root execution bucket for Cromwell
    */
  private def checkGcs(): Future[SubsystemStatus] = {
    val gcsAuthName = papiConfig.as[String]("filesystems.gcs.auth")
    val rootExecutionBucket = papiConfig.as[String]("root").replace("gs://", "")

    val cred = googleConfig.auth(gcsAuthName) match {
      case Valid(a) => a.credential
      case Invalid(e) => throw new IllegalArgumentException("Unable to configure WorkbenchHealthMonitor: " + e.toList.mkString(", "))
    }

    val storage = cred map { c => GcsStorage.gcsStorage(googleConfig.applicationName, c) }
    storage map { _.buckets.get(rootExecutionBucket).execute() } map { _ => OkStatus }
  }

  /**
    * Demonstrates connectivity to Google Pipelines API (PAPI) by making sure it can access an authenticated endpoint
    */
  private def checkPapi(): Future[SubsystemStatus] = {
    val endpointUrl = new URL(papiConfig.as[String]("genomics.endpoint-url"))
    val papiProjectId = papiConfig.as[String]("project")
    val papiAuthName = papiConfig.as[String]("genomics.auth")
    val papiAuthMode = googleConfig.auth(papiAuthName) match {
      case Valid(a) => a
      case Invalid(e) => throw new IllegalArgumentException("Unable to configure WorkbenchHealthMonitor: " + e.toList.mkString(", "))
    }

    val genomicsFactory = GenomicsFactory(googleConfig.applicationName, papiAuthMode, endpointUrl)
    val genomicsInterface = papiAuthMode.credential() map genomicsFactory.fromCredentials

    genomicsInterface map { gi =>
      gi.pipelines().list().setProjectId(papiProjectId).setPageSize(1).execute()
    } map { _ => OkStatus }
  }
}
