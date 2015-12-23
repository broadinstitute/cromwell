package cromwell.webservice

import cromwell.CromwellTestkitSpec
import cromwell.webservice.WorkflowQueryKey._
import org.joda.time.DateTime
import scalaz.{Name => _, _}

class WorkflowQueryParametersSpec extends CromwellTestkitSpec {

  val StartDateString = "2015-11-01T11:11:11"
  val EndDateString = "2015-11-01T12:12:12"

  "Workflow query parameters" should {

    "be accepted if empty" in {
      val result = WorkflowQueryParameters.runValidation(Seq.empty)
      result match {
        case Success(r) =>
          r.startDate should be('empty)
          r.endDate should be('empty)
          r.names should be('empty)
          r.statuses should be('empty)
        case Failure(fs) =>
          throw new RuntimeException(fs.list.mkString(", "))
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
        case Success(r) =>
          r.startDate.get should equal(new DateTime(StartDateString))
          r.endDate.get should equal(new DateTime(EndDateString))
          r.names should be(Set("my_workflow", "my_other_workflow"))
          r.statuses should be(Set("Succeeded", "Running"))
        case Failure(fs) =>
          throw new RuntimeException(fs.list.mkString(", "))
      }
    }

    "reject unrecognized keys" in {
      val rawParameters = Seq(
        Name.name -> "my_little_workflow",
        "Bogosity" -> ""
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Success(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Failure(fs) =>
          fs.list should have size 1
          fs.list.head should include("Unrecognized query keys: Bogosity")
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
        case Success(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Failure(fs) =>
          fs.list should have size 1
          fs.list.head should include("Specified start date is after specified end date")
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
        case Success(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Failure(fs) =>
          fs.list should have size 1
          fs.list.head should include("Name values do not match allowed workflow naming pattern")
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
        case Success(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Failure(fs) =>
          fs.list should have size 1
          fs.list.head should include("does not parse as a datetime")
      }
    }

    "reject multiple dates" in {
      val rawParameters = Seq(
        StartDate.name -> StartDateString,
        StartDate.name -> StartDateString
      )
      val result = WorkflowQueryParameters.runValidation(rawParameters)
      result match {
        case Success(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Failure(fs) =>
          fs.list should have size 1
          fs.list.head should include("at most one is allowed")
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
        case Success(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Failure(fs) =>
          fs.list should have size 1
          fs.list.head should be("Unrecognized status values: Moseying")
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
        case Success(r) =>
          throw new RuntimeException(s"Unexpected success: $r")
        case Failure(fs) =>
          fs.list should have size 3
          fs.list find { _ == "Unrecognized status values: Moseying" } getOrElse fail
          fs.list find { _ contains "does not parse as a datetime" } getOrElse fail
          fs.list find { _ contains "Name values do not match allowed workflow naming pattern" } getOrElse fail
      }
    }
  }
}
