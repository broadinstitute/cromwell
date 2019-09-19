package cromwell.services.database

import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

class DeleteMetadataRowsSpec extends FlatSpec with Matchers with ScalaFutures {

  // attempt deleting subworkflow id. should error

  // attempt deleting a subworkflow that has subworkflows itself. should error

  // attempt deleting root workflow id w/o subworkflows. should delete expected row count

  // attempt deleting root workflow id with subworkflows. should delete expected row count

  // attempt deleting root workflow id with subworkflows that have subworkflows. should delete expected row count

}
