package cromwell.services.healthmonitor.impl.workbench

import java.net.URL

import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.Config
import cromwell.cloudSupport.gcp.GoogleConfiguration
import cromwell.cloudSupport.gcp.gcs.GcsStorage
import cromwell.cloudSupport.gcp.genomics.GenomicsFactory
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
  val papiBackendName = serviceConfig.as[Option[String]]("services.HealthMonitor.PapiBackendName").getOrElse("Jes")
  val papiConfig: Config = globalConfig.as[Config]("backend.providers.Jes.config")

  /**
    * Demonstrates connectivity to GCS by stat-ing a public bucket. We use a public bucket as our plethora of auth
    * modes make it impossible to ensure we would always be able to ping an arbitrary bucket with the credentials
    * available.
    */
  private def checkGcs(): Future[SubsystemStatus] = {
    import WorkbenchHealthMonitorServiceActor.GcsBucketToCheck
    val gcsAuthName = papiConfig.as[String]("filesystems.gcs.auth")

    val cred = googleConfig.auth(gcsAuthName) match {
      case Valid(a) => a.credential(WorkflowOptions.empty)
      case Invalid(e) => throw new IllegalArgumentException("Unable to configure WorkbenchHealthMonitor: " + e.toList.mkString(", "))
    }

    val storage = cred map { c => GcsStorage.gcsStorage(googleConfig.applicationName, c) }
    storage map { _.buckets.get(GcsBucketToCheck).execute() } as OkStatus
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
    val genomicsInterface = papiAuthMode.credential(WorkflowOptions.empty) map genomicsFactory.fromCredentials

    genomicsInterface map { _.pipelines().list().setProjectId(papiProjectId).setPageSize(1).execute() } as OkStatus
  }
}

object WorkbenchHealthMonitorServiceActor {
  /*
    A bucket which is publicly readable. We might be authing to GCS with per-user credentials which means that we don't
    have generic Cromwell-level access to any private well known location.
   */
  val GcsBucketToCheck = "cromwell-health-check"
}
