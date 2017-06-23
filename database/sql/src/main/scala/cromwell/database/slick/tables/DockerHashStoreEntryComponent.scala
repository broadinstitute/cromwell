package cromwell.database.slick.tables

import cromwell.database.sql.tables.DockerHashStoreEntry

trait DockerHashStoreEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class DockerHashStoreEntries(tag: Tag) extends Table[DockerHashStoreEntry](tag, "DOCKER_HASH_STORE_ENTRY") {
    def dockerHashStoreEntryId = column[Int]("DOCKER_HASH_STORE_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def dockerTag = column[String]("DOCKER_TAG")

    def dockerHash = column[String]("DOCKER_HASH")

    override def * = (workflowExecutionUuid, dockerTag, dockerHash, dockerHashStoreEntryId.?) <> (DockerHashStoreEntry.tupled, DockerHashStoreEntry.unapply)

    def ucDockerHashStoreEntryWeuDt = index("UC_DOCKER_HASH_STORE_ENTRY_WEU_DT", (workflowExecutionUuid, dockerTag), unique = true)
  }

  val dockerHashStoreEntries = TableQuery[DockerHashStoreEntries]

  val dockerHashStoreEntryIdsAutoInc = dockerHashStoreEntries returning dockerHashStoreEntries.map(_.dockerHashStoreEntryId)

  /**
    * Useful for finding the docker hash store for a given workflow execution UUID
    */
  val dockerHashStoreEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      dockerHashStoreEntry <- dockerHashStoreEntries
      if dockerHashStoreEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield dockerHashStoreEntry
  )
}
