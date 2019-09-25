package cromwell.services.database

import cromwell.core.Tags.DbmsTest
import org.scalatest.concurrent.PatienceConfiguration.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.concurrent.ExecutionContext

class MetadataSlickDatabaseSpec extends FlatSpec with Matchers with ScalaFutures {

  DatabaseSystem.All foreach { databaseSystem =>

    implicit val ec = ExecutionContext.global

    behavior of s"MetadataSlickDatabase on ${databaseSystem.name}"

    lazy val dataAccess = DatabaseTestKit.initializedDatabaseFromSystem(MetadataDatabaseType, databaseSystem)

    // attempt deleting root workflow that does not exist
    it should "error when deleting a root workflow that does not exist" taggedAs DbmsTest in {
      val delete = dataAccess.deleteNonLabelMetadataForWorkflow("asdf")

      delete.failed.futureValue(Timeout(10.seconds)).getMessage should be(s"Programmer error: attempted to delete metadata rows for non-existent root workflow asdf")
    }

    // attempt deleting subworkflow id. should error

    // attempt deleting a subworkflow that has subworkflows itself. should error

    // attempt deleting root workflow id w/o subworkflows. should delete expected row count

    // attempt deleting root workflow id with subworkflows. should delete expected row count

    // attempt deleting root workflow id with subworkflows that have subworkflows. should delete expected row count

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.close()
    }
  }


}
