package cromwell.services.healthmonitor.impl.workbench

import java.net.URL

import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.Config
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.cloudsupport.gcp.genomics.GenomicsFactory
import cromwell.core.WorkflowOptions
import cromwell.services.healthmonitor.HealthMonitorServiceActor
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{MonitoredSubsystem, OkStatus, SubsystemStatus}
import net.ceedubs.ficus.Ficus._
import cats.instances.future._
import cats.syntax.functor._
import cromwell.services.healthmonitor.impl.common.{DockerHubMonitor, EngineDatabaseMonitor}

import scala.concurrent.Future

/**
  * A health monitor implementation which will monitor external services used in the Workbench ecosystem, such
  * as GCS and PAPI. This implementation makes some assumptions of Cromwell's configuration which will be true
  * in a Workbench scenario but YMMV otherwise. Caveat emptor and all of that fun stuff.
  */
class WorkbenchHealthMonitorServiceActor(val serviceConfig: Config, globalConfig: Config)
  extends HealthMonitorServiceActor
    with DockerHubMonitor
    with EngineDatabaseMonitor {
  override implicit val system = context.system

  override lazy val subsystems: Set[MonitoredSubsystem] = Set(DockerHub, EngineDb, Papi, Gcs)

  private lazy val Gcs = MonitoredSubsystem("GCS", checkGcs _)
  private lazy val Papi = MonitoredSubsystem("PAPI", checkPapi _)

  val googleConfig = GoogleConfiguration(globalConfig)

  val papiBackendName = serviceConfig.as[Option[String]]("papi-backend-name").getOrElse("Jes")
  val papiConfig: Config = globalConfig.as[Config](s"backend.providers.$papiBackendName.config")

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
    val storage = googleAuth.credential(WorkflowOptions.empty) map { c => GcsStorage.gcsStorage(googleConfig.applicationName, c) }
    storage map { _.buckets.get(gcsBucketToCheck).execute() } as OkStatus
  }

  /**
    * Demonstrates connectivity to Google Pipelines API (PAPI) by making sure it can access an authenticated endpoint
    */
  private def checkPapi(): Future[SubsystemStatus] = {
    val endpointUrl = new URL(papiConfig.as[String]("genomics.endpoint-url"))
    val papiProjectId = papiConfig.as[String]("project")

    val genomicsFactory = GenomicsFactory(googleConfig.applicationName, googleAuth, endpointUrl)
    val genomicsInterface = googleAuth.credential(WorkflowOptions.empty) map genomicsFactory.fromCredentials

    genomicsInterface map { _.pipelines().list().setProjectId(papiProjectId).setPageSize(1).execute() } as OkStatus
  }
}

