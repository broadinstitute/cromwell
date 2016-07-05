package cromwell.jobstore

import cromwell.core._

sealed trait JobResult

case class JobResultSuccess(returnCode: Option[Int], jobOutputs: JobOutputs) extends JobResult
case class JobResultFailure(returnCode: Option[Int], reason: Throwable) extends JobResult
