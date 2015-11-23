package cromwell.webservice

import cromwell.engine.WorkflowState
import org.joda.time.DateTime

import scala.language.{postfixOps, reflectiveCalls}
import scala.util.{Success, Try}
import scalaz.Scalaz._
import scalaz.{Success => _, _}

object WorkflowQueryKey {
  val ValidKeys = Set(StartDate, EndDate, Name, Status) map { _.name }

  case object StartDate extends DateTimeWorkflowQueryKey {
    override val name = "Start"
    override def displayName = "start date"
  }

  case object EndDate extends DateTimeWorkflowQueryKey {
    override val name = "End"
    override def displayName = "end date"
  }

  case object Name extends SeqStringWorkflowQueryKey {
    override val name = "Name"
    private val WorkflowNamePattern = "([a-zA-Z][a-zA-Z0-9_]+)".r

    override def validate(grouped: Map[String, Seq[(String, String)]]): ValidationNel[String, Seq[String]] = {
      val values = valuesFromMap(grouped).toList
      val nels = values map {
        case WorkflowNamePattern(n) => n.successNel
        case v => v.failureNel
      }
      sequenceListOfValidationNels("Name values do not match allowed workflow naming pattern [a-zA-Z][a-zA-Z0-9_]+", nels)
    }
  }

  case object Status extends SeqStringWorkflowQueryKey {
    override val name = "Status"

    override def validate(grouped: Map[String, Seq[(String, String)]]): ValidationNel[String, Seq[String]] = {
      val values = valuesFromMap(grouped).toList
      val nels = values map { v =>
        if (Try(WorkflowState.fromString(v.toLowerCase.capitalize)).isSuccess) v.successNel else v.failureNel
      }
      sequenceListOfValidationNels("Unrecognized status values", nels)
    }
  }
}

sealed trait WorkflowQueryKey[T] {
  def validate(grouped: Map[String, Seq[(String, String)]]): ValidationNel[String, T]
  def name: String
  def valuesFromMap(grouped: Map[String, Seq[(String, String)]]): Seq[String] = {
    grouped.getOrElse(name, Seq.empty) map { _._2 }
  }
}

sealed trait DateTimeWorkflowQueryKey extends WorkflowQueryKey[Option[DateTime]] {
  override def validate(grouped: Map[String, Seq[(String, String)]]): ValidationNel[String, Option[DateTime]] = {
    valuesFromMap(grouped) match {
      case vs if vs.size > 1 =>
        s"Found ${vs.size} values for key '$name' but at most one is allowed.".failureNel
      case Nil => None.successNel
      case v :: Nil =>
        Try(new DateTime(v)) match {
          case Success(dt) => Option(dt).successNel
          case _ => s"Value given for $displayName does not parse as a datetime: $v".failureNel
        }
    }
  }
  def displayName: String
}

sealed trait SeqStringWorkflowQueryKey extends WorkflowQueryKey[Seq[String]] {

  /** `sequence` the `List[ValidationNel[String, String]]` to a single `ValidationNel[String, List[String]]` */
  protected def sequenceListOfValidationNels(prefix: String,
                                             validationNelList: List[ValidationNel[String, String]]): ValidationNel[String, List[String]] = {

    val validationNel = validationNelList.sequence[({type λ[a] = ValidationNel[String, a]})#λ, String]
    // With a leftMap, prepend an error message to the concatenated error values if there are error values.
    // This turns the ValidationNel into a Validation, force it back to a ValidationNel with toValidationNel.
    validationNel.leftMap(prefix + ": " + _.list.mkString(", ")).toValidationNel
  }
}

