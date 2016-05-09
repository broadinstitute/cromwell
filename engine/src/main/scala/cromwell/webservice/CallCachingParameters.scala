package cromwell.webservice

import cromwell.core.{ErrorOr, WorkflowId}
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.{DataAccess, ExecutionDatabaseKey}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try
import scalaz.Scalaz._
import scalaz.{Failure, Success}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CallCachingParameters private(workflowId: WorkflowId, callKey: Option[ExecutionDatabaseKey], allow: Boolean)

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object CallCachingParameters {
  private [webservice] def validateRecognizedKeys(queryParameters: QueryParameters): ErrorOr[Unit] = {
    val badKeys = queryParameters collect { case q if q.key.toLowerCase != "allow" => q.key }
    if (badKeys.nonEmpty) ("Found unrecognized keys: " + badKeys.mkString(", ")).failureNel else ().successNel
  }

  private [webservice] def validateAllow(queryParameters: QueryParameters): ErrorOr[Boolean] = {
    val allows = queryParameters.collect({ case q if q.key.toLowerCase == "allow" => q.value }).toSet
    val (booleans, nonBooleans) = allows map { a => (a, Try(a.toBoolean)) } partition { _._2.isSuccess }

    val allBooleans = if (nonBooleans.nonEmpty) {
      val values = (nonBooleans map { _._1 } ).mkString(", ")
      s"Found non-boolean 'allow' values: $values".failureNel
    } else {
      ().successNel
    }

    val (trues, falses) = booleans partition { _._2.get }
    val coherentValues = if (trues.nonEmpty && falses.nonEmpty) "Found both true and false 'allow' values".failureNel else trues.nonEmpty.successNel

    (allBooleans |@| coherentValues) {
      case (_, allow) => allow
    }
  }

  private def validateWorkflowExists(workflowId: WorkflowId, dataAccess: DataAccess)(implicit ec: ExecutionContext): Future[ErrorOr[Unit]] = {
    dataAccess.getWorkflowState(workflowId) map {
      case Some(_) => ().successNel
      case None => s"Workflow not found: ${workflowId.id.toString}".failureNel
    }
  }

  private [webservice] def validateCallName(callName: Option[String]): ErrorOr[Option[ExecutionDatabaseKey]] = {
    import Patterns.CallFullyQualifiedName
    callName map {
      case CallFullyQualifiedName(fqn, index) => Option(ExecutionDatabaseKey(fqn, Option(index) map { _.toInt }, 1)).successNel
      case name => s"Specified call does not parse as a fully qualified name with optional index: $name".failureNel
    } getOrElse None.successNel
  }

  /**
   * Assert that a call with the specified key exists in the workflow.
   */
  private def validateCall(workflowId: WorkflowId, key: ExecutionDatabaseKey, dataAccess: DataAccess)
                          (implicit ec: ExecutionContext): Future[ErrorOr[Unit]] = {

    dataAccess.getExecutionStatus(workflowId, key) map {
      case Some(s) => ().successNel
      case None =>
        val displayIndex = key.index map { "." + _ } getOrElse ""
        s"Call ${key.fqn}$displayIndex does not exist in workflow $workflowId".failureNel
    }
  }

  private def validateWorkflowAndCall(workflowId: WorkflowId, callName: Option[String], dataAccess: DataAccess)
                                     (implicit ec: ExecutionContext): Future[ErrorOr[Option[ExecutionDatabaseKey]]] = {
    /**
     * Perform a call validation conditioned on the workflow and call name validations having already succeeded.
     * There's no possibility that the particulars of the call are even in the database if the workflow ID doesn't
     * exist or the call name is malformed.
     */
    def validateCallConditionally(workflowValidation: ErrorOr[Unit],
                                  callNameValidation: ErrorOr[Option[ExecutionDatabaseKey]],
                                  dataAccess: DataAccess)
                                 (implicit ec: ExecutionContext): Future[ErrorOr[Unit]] = {

      (workflowValidation, callNameValidation) match {
        case (Success(_), Success(call)) =>
          call map { c => validateCall(workflowId, c, dataAccess) } getOrElse Future.successful(().successNel)
        // If either the workflow or call name validations are failed don't produce additional validation failures.
        case _ => Future.successful(().successNel)
      }
    }

    for {
      workflowValidation <- validateWorkflowExists(workflowId, dataAccess)
      callNameValidation <- Future(validateCallName(callName))
      callValidation <- validateCallConditionally(workflowValidation, callNameValidation, dataAccess)
    } yield (workflowValidation |@| callNameValidation |@| callValidation) {
      case (_, key, _) => key
    }
  }

  def from(workflowId: WorkflowId, callName: Option[String], queryParameters: QueryParameters, dataAccess: DataAccess = globalDataAccess)
          (implicit ec: ExecutionContext): Future[CallCachingParameters] = {

    val validations = for {
      recognizedKeys <- Future(validateRecognizedKeys(queryParameters))
      validWorkflowAndCall <- validateWorkflowAndCall(workflowId, callName, dataAccess)
      allow <- Future(validateAllow(queryParameters))
    } yield (recognizedKeys |@| validWorkflowAndCall |@| allow) {
      case (_, callKey, a) => CallCachingParameters(workflowId, callKey, a)
    }

    validations map {
      case Success(s) => s
      case Failure(errors) => throw new IllegalArgumentException(errors.list.mkString("\n"))
    }
  }
}
