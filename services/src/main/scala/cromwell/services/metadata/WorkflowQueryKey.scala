package cromwell.services.metadata

import java.time.OffsetDateTime

import cats.data
import cats.syntax.traverse._
import cats.syntax.validated._
import cromwell.core.labels.Label
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState}
import common.validation.ErrorOr._
import cats.data.Validated._
import cats.instances.list._
import mouse.boolean._

import scala.util.{Success, Try}

object WorkflowQueryKey {
  val ValidKeys = Set(
    StartDate,
    EndDate,
    Name,
    Id,
    Status,
    LabelAndKeyValue,
    LabelOrKeyValue,
    ExcludeLabelAndKeyValue,
    ExcludeLabelOrKeyValue,
    Page,
    PageSize,
    AdditionalQueryResultFields,
    SubmissionTime,
    IncludeSubworkflows
  ) map { _.name }

  case object StartDate extends DateTimeWorkflowQueryKey {
    override val name = "Start"
    override def displayName = "start date"
  }

  case object EndDate extends DateTimeWorkflowQueryKey {
    override val name = "End"
    override def displayName = "end date"
  }

  case object SubmissionTime extends DateTimeWorkflowQueryKey {
    override val name = "Submission"
    override def displayName = "submission time"
  }

  case object Page extends IntWorkflowQueryKey {
    override val name = "Page"
    override def displayName = "page"
  }

  case object PageSize extends IntWorkflowQueryKey {
    override val name = "Pagesize"
    override def displayName = "page size"
  }

  case object Name extends SeqWorkflowQueryKey[String] {
    override val name = "Name"

    override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[List[String]] = {
      val values = valuesFromMap(grouped).toList
      val nels:List[data.ValidatedNel[String,String]] = values map {
        case Patterns.WorkflowName(n) => n.validNel[String]
        case v => v.invalidNel[String]
      }
      sequenceListOfValidatedNels("Name values do not match allowed workflow naming pattern", nels)
    }
  }

  sealed trait LabelLikeKeyValue extends SeqWorkflowQueryKey[Label] {
    override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[List[Label]] = {
      val values = valuesFromMap(grouped).toList

      def validateLabelRegex(labelKeyValue: String): ErrorOr[Label] = {
        labelKeyValue.split("\\:", 2) match {
          case Array(k, v) => Label.validateLabel(k, v)
          case _ => labelKeyValue.invalidNel
        }
      }
      val nels: List[ErrorOr[Label]] = values map validateLabelRegex
      sequenceListOfValidatedNels("Label values do not match allowed pattern label-key:label-value", nels)
    }
  }


  case object LabelAndKeyValue extends LabelLikeKeyValue {
    override val name = "Label"
  }

  case object LabelOrKeyValue extends LabelLikeKeyValue {
    override val name = "Labelor"
  }

  case object ExcludeLabelAndKeyValue extends LabelLikeKeyValue {
    override val name = "Excludelabeland"
  }

  case object ExcludeLabelOrKeyValue extends LabelLikeKeyValue {
    override val name = "Excludelabelor"
  }

  case object Id extends SeqWorkflowQueryKey[String] {
    override val name = "Id"

    override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[List[String]] = {
      val values = valuesFromMap(grouped).toList
      val nels = values map { v =>
        if (Try(WorkflowId.fromString(v.toLowerCase.capitalize)).isSuccess) v.validNel[String] else v.invalidNel[String]
      }
      sequenceListOfValidatedNels("Id values do match allowed workflow id pattern", nels)
    }
  }

  case object Status extends SeqWorkflowQueryKey[String] {
    override val name = "Status"

    override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[List[String]] = {
      val values = valuesFromMap(grouped).toList
      val nels = values map { v =>
        Try(WorkflowState.withName(v)) match {
          case Success(workflowState) => workflowState.toString.validNel[String]
          case _ => v.invalidNel[String]
        }
      }
      sequenceListOfValidatedNels("Unrecognized status values", nels)
    }
  }

  case object AdditionalQueryResultFields extends SeqWorkflowQueryKey[String] {
    override val name = "Additionalqueryresultfields"

    override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[List[String]] = {
      val values = valuesFromMap(grouped).toList
      /*
        The inclusion of `WorkflowMetadataKeys.ParentWorkflowId` is for backwards compatibility. As of #4381
        parentWorkflowId is always included, but we did not want to break old automated queries
        */
      val allowedValues = Seq(WorkflowMetadataKeys.Labels, WorkflowMetadataKeys.ParentWorkflowId)
      val nels: List[ErrorOr[String]] = values map { v => {
        allowedValues.contains(v).fold(v.validNel[String], v.invalidNel[String])
      }}
      sequenceListOfValidatedNels(s"Keys should be from $allowedValues. Unrecognized values", nels)
    }
  }

  case object IncludeSubworkflows extends BooleanWorkflowQueryKey {
    override val name = "Includesubworkflows"
    override def displayName = "include subworkflows"
    override def defaultBooleanValue: Boolean = true
  }
}

sealed trait WorkflowQueryKey[T] {
  def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[T]
  def name: String
  def valuesFromMap(grouped: Map[String, Seq[(String, String)]]): Seq[String] = {
    grouped.getOrElse(name, Seq.empty) map { _._2 }
  }
}

sealed trait DateTimeWorkflowQueryKey extends WorkflowQueryKey[Option[OffsetDateTime]] {
  override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[Option[OffsetDateTime]] = {
    valuesFromMap(grouped).toList match {
      case vs if vs.lengthCompare(1) > 0 =>
        s"Found ${vs.size} values for key '$name' but at most one is allowed.".invalidNel[Option[OffsetDateTime]]
      case Nil => None.validNel[String]
      case v :: Nil =>
        Try(OffsetDateTime.parse(v)) match {
          case Success(dt) => Option(dt).validNel[String]
          case _ => s"Value given for $displayName does not parse as a datetime: $v".invalidNel[Option[OffsetDateTime]]
        }
    }
  }
  def displayName: String
}

sealed trait SeqWorkflowQueryKey[A] extends WorkflowQueryKey[Seq[A]] {
  /** `sequence` the `List[ErrorOr[A]]` to a single `ErrorOr[List[A]]` */
  protected def sequenceListOfValidatedNels(prefix: String, errorOrList: List[ErrorOr[A]]): ErrorOr[List[A]] = {
    val errorOr = errorOrList.sequence[ErrorOr, A]
    // With a leftMap, prepend an error message to the concatenated error values if there are error values.
    // This turns the ValidatedNel into a Validated, force it back to a ValidatedNel with toValidationNel.
    errorOr.leftMap(prefix + ": " + _.toList.mkString(", ")).toValidatedNel
  }
}

sealed trait IntWorkflowQueryKey extends WorkflowQueryKey[Option[Int]] {
  override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[Option[Int]] = {
    valuesFromMap(grouped).toList match {
      case vs if vs.lengthCompare(1) > 0 =>
        s"Found ${vs.size} values for key '$name' but at most one is allowed.".invalidNel[Option[Int]]
      case Nil => None.validNel
      case v :: Nil =>
        Try(v.toInt) match {
          case Success(intVal) => if (intVal > 0) Option(intVal).validNel else s"Integer value not greater than 0".invalidNel[Option[Int]]
          case _ => s"Value given for $displayName does not parse as a integer: $v".invalidNel[Option[Int]]
        }
    }
  }
  def displayName: String
}

sealed trait BooleanWorkflowQueryKey extends WorkflowQueryKey[Boolean] {
  override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[Boolean] = {
    valuesFromMap(grouped).toList match {
      case vs if vs.lengthCompare(1) > 0 => s"Found ${vs.size} values for key '$name' but at most one is allowed.".invalidNel[Boolean]
      case Nil => defaultBooleanValue.validNel
      case v :: Nil => {
        Try(v.toBoolean) match {
          case Success(bool) => bool.validNel
          case _ => s"Value given for $displayName does not parse as a boolean: $v".invalidNel[Boolean]
        }
      }
    }
  }
   def displayName: String

  def defaultBooleanValue: Boolean
}

