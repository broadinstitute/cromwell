package cromwell.services.metadata

import java.time.OffsetDateTime

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import cromwell.core.{WorkflowId, WorkflowState}
import lenthall.validation.ErrorOr._

import scala.util.{Success, Try}

object WorkflowQueryKey {
  val ValidKeys = Set(StartDate, EndDate, Name, Id, Status, Page, PageSize) map { _.name }

  case object StartDate extends DateTimeWorkflowQueryKey {
    override val name = "Start"
    override def displayName = "start date"
  }

  case object EndDate extends DateTimeWorkflowQueryKey {
    override val name = "End"
    override def displayName = "end date"
  }

  case object Page extends IntWorkflowQueryKey {
    override val name = "Page"
    override def displayName = "page"
  }

  case object PageSize extends IntWorkflowQueryKey {
    override val name = "Pagesize"
    override def displayName = "page size"
  }

  case object Name extends SeqStringWorkflowQueryKey {
    override val name = "Name"

    override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[List[String]] = {

      val values = valuesFromMap(grouped).toList
      val nels = values map {
        case Patterns.WorkflowName(n) => n.validNel[String]
        case v => v.invalidNel[String]
      }
      sequenceListOfValidatedNels(s"Name values do not match allowed workflow naming pattern", nels)
    }
  }

  case object Id extends SeqStringWorkflowQueryKey {
    override val name = "Id"

    override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[List[String]] = {
      val values = valuesFromMap(grouped).toList
      val nels = values map { v =>
        if (Try(WorkflowId.fromString(v.toLowerCase.capitalize)).isSuccess) v.validNel[String] else v.invalidNel[String]
      }
      sequenceListOfValidatedNels(s"Id values do match allowed workflow id pattern", nels)
    }
  }

  case object Status extends SeqStringWorkflowQueryKey {
    override val name = "Status"

    override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[List[String]] = {
      val values = valuesFromMap(grouped).toList
      val nels = values map { v =>
        if (Try(WorkflowState.fromString(v.toLowerCase.capitalize)).isSuccess) v.validNel[String] else v.invalidNel[String]
      }
      sequenceListOfValidatedNels("Unrecognized status values", nels)
    }
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
    valuesFromMap(grouped) match {
      case vs if vs.size > 1 =>
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

sealed trait SeqStringWorkflowQueryKey extends WorkflowQueryKey[Seq[String]] {
  /** `sequence` the `List[ErrorOr[String]]` to a single `ErrorOr[List[String]]` */
  protected def sequenceListOfValidatedNels(prefix: String, errorOrList: List[ErrorOr[String]]): ErrorOr[List[String]] = {
    val errorOr = errorOrList.sequence[ErrorOr, String]
    // With a leftMap, prepend an error message to the concatenated error values if there are error values.
    // This turns the ValidatedNel into a Validated, force it back to a ValidatedNel with toValidationNel.
    errorOr.leftMap(prefix + ": " + _.toList.mkString(", ")).toValidatedNel
  }
}

sealed trait IntWorkflowQueryKey extends WorkflowQueryKey[Option[Int]] {
  override def validate(grouped: Map[String, Seq[(String, String)]]): ErrorOr[Option[Int]] = {
    valuesFromMap(grouped) match {
      case vs if vs.size > 1 =>
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

