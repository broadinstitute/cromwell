package centaur.reporting

import java.time.OffsetDateTime
import java.util

import cats.effect.IO
import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.traverse._
import centaur.reporting.BigQueryReporter._
import centaur.test.CentaurTestException
import centaur.test.metadata.CallAttemptFailure
import com.google.api.gax.retrying.RetrySettings
import com.google.api.services.bigquery.BigqueryScopes
import com.google.auth.Credentials
import com.google.cloud.bigquery.InsertAllRequest.RowToInsert
import com.google.cloud.bigquery.{BigQuery, BigQueryError, BigQueryOptions, InsertAllRequest, InsertAllResponse, TableId}
import common.util.TimeUtil._
import common.validation.Validation._
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.tables.{JobKeyValueEntry, MetadataEntry}
import net.ceedubs.ficus.Ficus._
import org.apache.commons.lang3.exception.ExceptionUtils
import org.threeten.bp.Duration

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext

class BigQueryReporter(override val params: ErrorReporterParams) extends ErrorReporter {

  private val bigQueryAuth: String = params.reporterConfig.getString("auth")
  private val bigQueryProjectOption: Option[String] = params.reporterConfig.getAs[String]("project")
  private val bigQueryDataset: String = params.reporterConfig.getString("dataset")

  override lazy val destination: String = bigQueryProjectOption.map(_ + "/").getOrElse("") + bigQueryDataset

  private val retrySettings: RetrySettings = {
    RetrySettings.newBuilder()
      .setMaxAttempts(3)
      .setTotalTimeout(Duration.ofSeconds(30))
      .setInitialRetryDelay(Duration.ofMillis(100))
      .setRetryDelayMultiplier(1.1)
      .setMaxRetryDelay(Duration.ofSeconds(1))
      .setInitialRpcTimeout(Duration.ofMillis(100))
      .setRpcTimeoutMultiplier(1.1)
      .setMaxRpcTimeout(Duration.ofSeconds(5))
      .build()
  }

  private val bigQueryCredentials: Credentials = GoogleConfiguration
    .apply(params.rootConfig)
    .auth(bigQueryAuth)
    .unsafe
    .credentials(Set(BigqueryScopes.BIGQUERY_INSERTDATA))

  private val bigQuery: BigQuery = BigQueryOptions.newBuilder()
    .setRetrySettings(retrySettings)
    .setCredentials(bigQueryCredentials)
    .build()
    .getService

  private val testFailureTableId = bigQueryTable("test_failure")
  private val callAttemptFailureTableId = bigQueryTable("call_attempt_failure")
  private val jobKeyValueTableId = bigQueryTable("job_key_value")
  private val metadataTableId = bigQueryTable("metadata")

  def bigQueryTable(table: String): TableId = {
    bigQueryProjectOption match {
      case Some(project) => TableId.of(project, bigQueryDataset, table)
      case None => TableId.of(bigQueryDataset, table)
    }
  }

  override def logCentaurFailure(testEnvironment: TestEnvironment,
                                 ciEnvironment: CiEnvironment,
                                 centaurTestException: CentaurTestException)
                                (implicit executionContext: ExecutionContext): IO[Unit] = {
    for {
      callAttemptFailures <- CallAttemptFailure.buildFailures(centaurTestException.metadataJsonOption)
      jobKeyValueEntries <- params.database.jobKeyValueEntriesIo(centaurTestException.workflowIdOption)
      metadataEntries <- params.database.metadataEntriesIo(centaurTestException.workflowIdOption)
      _ <- sendBigQueryFailure(
        testEnvironment,
        ciEnvironment,
        centaurTestException,
        callAttemptFailures,
        jobKeyValueEntries,
        metadataEntries
      )
    } yield ()
  }

  private def sendBigQueryFailure(testEnvironment: TestEnvironment,
                                  ciEnvironment: CiEnvironment,
                                  centaurTestException: CentaurTestException,
                                  callAttemptFailures: Vector[CallAttemptFailure],
                                  jobKeyValueEntries: Seq[JobKeyValueEntry],
                                  metadataEntries: Seq[MetadataEntry]): IO[Unit] = {

    val metadata: IO[List[BigQueryError]] = {
      val metadataRows: List[util.List[RowToInsert]] = metadataEntries.map(toMetadataRow).grouped(10000).map(_.asJava).toList
      val metadataRequest: List[InsertAllRequest] = metadataRows.map(InsertAllRequest.of(metadataTableId, _))

      if (metadataEntries.nonEmpty)
        metadataRequest.traverse[IO, List[BigQueryError]](req => IO(bigQuery.insertAll(req)).map(_.getErrors)).map(_.flatten)
      else
        IO{Nil}
    }
    (IO {
      val testFailureRow = toTestFailureRow(testEnvironment, ciEnvironment, centaurTestException)
      val callAttemptFailureRows = callAttemptFailures.map(toCallAttemptFailureRow).asJava
      val jobKeyValueRows = jobKeyValueEntries.map(toJobKeyValueRow).asJava

      val testFailureRequest = InsertAllRequest.of(testFailureTableId, testFailureRow)
      val callAttemptFailuresRequest = InsertAllRequest.of(callAttemptFailureTableId, callAttemptFailureRows)
      val jobKeyValuesRequest = InsertAllRequest.of(jobKeyValueTableId, jobKeyValueRows)
      val testFailureErrors = bigQuery.insertAll(testFailureRequest).getErrors
      val callAttemptFailuresErrors =
        if (callAttemptFailures.nonEmpty) bigQuery.insertAll(callAttemptFailuresRequest).getErrors else Nil
      val jobKeyValuesErrors =
        if (jobKeyValueEntries.nonEmpty) bigQuery.insertAll(jobKeyValuesRequest).getErrors else Nil

      testFailureErrors ++ callAttemptFailuresErrors ++ jobKeyValuesErrors
    }, metadata).mapN(_ ++ _).flatMap {
      case errors if errors.isEmpty => IO.unit
      case errors => IO.raiseError {
        val errorCount = errors.size
        val threeErrors = errors.map(String.valueOf).distinct.sorted.take(3)
        val continued = if (errorCount > 3) "\n..." else ""
        val message = threeErrors.mkString(
          s"$errorCount error(s) occurred uploading to BigQuery: \n",
          "\n",
          continued)
        new RuntimeException(message)
      }
    }
  }

  private def toTestFailureRow(testEnvironment: TestEnvironment,
                               ciEnvironment: CiEnvironment,
                               centaurTestException: CentaurTestException): RowToInsert = {
    RowToInsert of Map(
      "ci_env_branch" -> ciEnvironment.branch,
      "ci_env_event" -> ciEnvironment.event,
      "ci_env_is_ci" -> ciEnvironment.isCi,
      "ci_env_number" -> ciEnvironment.number,
      "ci_env_os" -> ciEnvironment.os,
      "ci_env_provider" -> ciEnvironment.provider,
      "ci_env_tag" -> ciEnvironment.tag,
      "ci_env_type" -> ciEnvironment.`type`,
      "ci_env_url" -> ciEnvironment.url,
      "ci_env_centaur_type" -> ciEnvironment.centaurType,
      "test_attempt" -> Option(testEnvironment.attempt + 1),
      "test_message" -> Option(centaurTestException.message),
      "test_name" -> Option(testEnvironment.name),
      "test_stack_trace" -> Option(ExceptionUtils.getStackTrace(centaurTestException)),
      "test_timestamp" -> Option(OffsetDateTime.now.toUtcMilliString),
      "test_workflow_id" -> centaurTestException.workflowIdOption,
    ).collect {
      case (key, Some(value)) => (key, value)
    }.asJava
  }

  private def toCallAttemptFailureRow(callAttemptFailure: CallAttemptFailure): RowToInsert = {
    RowToInsert of Map(
      "call_fully_qualified_name" -> Option(callAttemptFailure.callFullyQualifiedName),
      "call_root" -> callAttemptFailure.callRootOption,
      "end" -> callAttemptFailure.endOption.map(_.toUtcMilliString),
      "job_attempt" -> Option(callAttemptFailure.jobAttempt),
      "job_index" -> Option(callAttemptFailure.jobIndex),
      "message" -> Option(callAttemptFailure.message),
      "start" -> callAttemptFailure.startOption.map(_.toUtcMilliString),
      "stderr" -> callAttemptFailure.stderrOption,
      "stdout" -> callAttemptFailure.stdoutOption,
      "workflow_id" -> Option(callAttemptFailure.workflowId),
    ).collect {
      case (key, Some(value)) => (key, value)
    }.asJava
  }

  private def toJobKeyValueRow(jobKeyValueEntry: JobKeyValueEntry): RowToInsert = {
    RowToInsert of Map(
      "call_fully_qualified_name" -> jobKeyValueEntry.callFullyQualifiedName,
      "job_attempt" -> jobKeyValueEntry.jobAttempt,
      "job_index" -> jobKeyValueEntry.jobIndex,
      "store_key" -> jobKeyValueEntry.storeKey,
      "store_value" -> jobKeyValueEntry.storeValue,
      "workflow_execution_uuid" -> jobKeyValueEntry.workflowExecutionUuid,
    ).asJava
  }

  private def toMetadataRow(metadataEntry: MetadataEntry): RowToInsert = {
    RowToInsert of Map(
      "call_fully_qualified_name" -> metadataEntry.callFullyQualifiedName,
      "job_attempt" -> metadataEntry.jobAttempt,
      "job_index" -> metadataEntry.jobIndex,
      "metadata_key" -> Option(metadataEntry.metadataKey),
      "metadata_timestamp" -> Option(metadataEntry.metadataTimestamp.toSystemOffsetDateTime.toUtcMilliString),
      "metadata_value" -> metadataEntry.metadataValue.map(_.toRawString),
      "metadata_value_type" -> metadataEntry.metadataValueType,
      "workflow_execution_uuid" -> Option(metadataEntry.workflowExecutionUuid),
    ).collect {
      case (key, Some(value)) => (key, value)
    }.asJava
  }
}

object BigQueryReporter {

  implicit class EnhancedInsertAllResponse(val response: InsertAllResponse) extends AnyVal {
    def getErrors: List[BigQueryError] = response.getInsertErrors.asScala.values.flatMap(_.asScala).toList
  }

}
