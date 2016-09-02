package cromwell.services.metadata

import java.time.OffsetDateTime

import cats.data.Validated._
import cromwell.services.metadata.WorkflowQueryKey._
import org.scalatest.{Matchers, WordSpec}

class WorkflowQueryParametersSpec extends WordSpec with Matchers {

  val StartDateString = "2015-11-01T11:11:11Z"
  val EndDateString = "2015-11-01T12:12:12Z"

  "Workflow query parameters" should {

    "be accepted if empty" in {
      val result = WorkflowQueryParameters.runValidation(Seq.empty)
      result match {
        case Valid(r) =>
          r.startDate should be('empty)
          r.endDate should be('empty)
          r.names should be('empty)
          r.statuses should be('empty)
        case Invalid(fs) =>
          throw new RuntimeException(fs.toList.mkString(", "))
      }
    }

    "accept correctly specified parameters" in {
      val rawParameters = Seq(
        Name.name -> "my_workflow",
        Status.name -> "Succeeded",
        Name.name -> "my_other_workflow",
        Status.name -> "Running",
        StartDate.name -> StartDateString,
        EndDate.name -> EndDateString
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Valid(r) =>
          r.startDate.get.toInstant should equal(OffsetDateTime.parse(StartDateString).toInstant)
          r.endDate.get.toInstant should equal(OffsetDateTime.parse(EndDateString).toInstant)
          r.names should be(Set("my_workflow", "my_other_workflow"))
          r.statuses should be(Set("Succeeded", "Running"))
        case Invalid(fs) =>
          throw new RuntimeException(fs.toList.mkString(", "))
      }
    }

    "reject unrecognized keys" in {
      val rawParameters = Seq(
        Name.name -> "my_little_workflow",
        "Bogosity" -> ""
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Valid(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Invalid(fs) =>
          fs.toList should have size 1
          fs.toList.head should include("Unrecognized query keys: Bogosity")
      }
    }

    "reject incompatible dates" in {
      // Intentionally setting start after end.
      val rawParameters = Seq(
        StartDate.name -> EndDateString,
        EndDate.name -> StartDateString
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Valid(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Invalid(fs) =>
          fs.toList should have size 1
          fs.toList.head should include("Specified start date is after specified end date")
      }
    }

    "reject nonconforming names" in {
      // One good one bad.
      val rawParameters = Seq(
        Name.name -> "my_workflow",
        // True statement but not at all suitable as an identifier.
        Name.name -> "WDL Rocks!!!"
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Valid(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Invalid(fs) =>
          fs.toList should have size 1
          fs.toList.head should include("Name values do not match allowed workflow naming pattern")
      }
    }

    "reject malformed dates" in {
      // One well-formed date and one badly-formed.
      val rawParameters = Seq(
        StartDate.name -> StartDateString,
        EndDate.name -> (EndDateString + "_")
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Valid(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Invalid(fs) =>
          fs.toList should have size 1
          fs.toList.head should include("does not parse as a datetime")
      }
    }

    "reject multiple dates" in {
      val rawParameters = Seq(
        StartDate.name -> StartDateString,
        StartDate.name -> StartDateString
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Valid(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Invalid(fs) =>
          fs.toList should have size 1
          fs.toList.head should include("at most one is allowed")
      }
    }

    "reject bad statuses" in {
      // One good one bad.
      val rawParameters = Seq(
        Status.name -> "Running",
        Status.name -> "Moseying"
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Valid(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Invalid(fs) =>
          fs.toList should have size 1
          fs.toList.head should be("Unrecognized status values: Moseying")
      }
    }

    "collect multiple errors for a hot mess of a query" in {
      val rawParameters = Seq(
        Status.name -> "Moseying",
        EndDate.name -> (EndDateString + "_"),
        Name.name -> "WDL Rocks!!!",
        StartDate.name -> StartDateString
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Valid(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Invalid(fs) =>
          fs.toList should have size 3
          fs.toList find { _ == "Unrecognized status values: Moseying" } getOrElse fail
          fs.toList find { _ contains "does not parse as a datetime" } getOrElse fail
          fs.toList find { _ contains "Name values do not match allowed workflow naming pattern" } getOrElse fail
      }
    }
  }
}
