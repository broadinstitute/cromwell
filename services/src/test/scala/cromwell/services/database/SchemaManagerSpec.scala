package cromwell.services.database

import java.io.{ByteArrayOutputStream, PrintStream}

import better.files._
import cromwell.database.migration.liquibase.LiquibaseUtils
import cromwell.database.slick.SlickDatabase
import cromwell.services.database.DatabaseTestKit._
import liquibase.diff.DiffResult
import liquibase.diff.output.DiffOutputControl
import liquibase.diff.output.changelog.DiffToChangeLog
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}
import slick.jdbc.JdbcProfile
import slick.jdbc.meta._

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.Try
import scala.xml._

/**
  * Tests that all table objects are consistently named using in-memory HSQLDB instances.
  */
class SchemaManagerSpec extends FlatSpec with Matchers with ScalaFutures {

  import SchemaManagerSpec._

  implicit val executionContext = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = scaled(5.seconds), interval = scaled(100.millis))

  private def getSchemaMetadata(slickDatabase: SlickDatabase): SchemaMetadata = {
    import slickDatabase.dataAccess.driver.api._
    val schemaMetadataFuture =
      for {
        tables <- slickDatabase.database.run(MTable.getTables(Option("PUBLIC"), Option("PUBLIC"), None, None))
        workingTables = tables
          .filterNot(_.name.name.contains("DATABASECHANGELOG"))
          // NOTE: MetadataEntry column names are perma-busted due to the large size of the table.
          .filterNot(_.name.name == "METADATA_ENTRY")
        columns <- slickDatabase.database.run(DBIO.sequence(workingTables.map(_.getColumns)))
        indexes <- slickDatabase.database.run(DBIO.sequence(workingTables.map(_.getIndexInfo())))
        primaryKeys <- slickDatabase.database.run(DBIO.sequence(workingTables.map(_.getPrimaryKeys)))
        foreignKeys <- slickDatabase.database.run(DBIO.sequence(workingTables.map(_.getExportedKeys)))
      } yield SchemaMetadata(
        tables,
        columns.flatten,
        indexes.flatten.filterNot(isGenerated),
        primaryKeys.flatten.filterNot(isGenerated),
        foreignKeys.flatten
      )

    schemaMetadataFuture.futureValue
  }

  CromwellDatabaseType.All foreach { databaseType =>
    SchemaManager.All foreach { schemaManager =>
      behavior of s"${databaseType.name} ${schemaManager.name}"

      val otherSchemaManager = schemaManager.other

      it should s"have the same schema as ${databaseType.name} ${otherSchemaManager.name}" in {
        for {
          actualDatabase <- inMemoryDatabase(databaseType, schemaManager).autoClosed
          expectedDatabase <- inMemoryDatabase(databaseType, otherSchemaManager).autoClosed
        } {
          compare(
            expectedDatabase.dataAccess.driver, expectedDatabase.database,
            actualDatabase.dataAccess.driver, actualDatabase.database,
          ) { diffResult =>

            import cromwell.database.migration.liquibase.DiffResultFilter._

            /*
            NOTE: Unique indexes no longer need to be filtered, as WE SHOULD NOT BE USING THEM!
            See notes at the bottom of changelog.xml
             */
            val diffFilters = StandardTypeFilters
            val filteredDiffResult = diffResult
              .filterChangeLogs
              .filterLiquibaseObjects
              .filterChangedObjects(diffFilters)

            val totalChanged =
              filteredDiffResult.getChangedObjects.size +
                filteredDiffResult.getMissingObjects.size +
                filteredDiffResult.getUnexpectedObjects.size

            if (totalChanged > 0) {
              val outputStream = new ByteArrayOutputStream
              val printStream = new PrintStream(outputStream, true)
              val diffOutputControl = new DiffOutputControl(false, false, false, Array.empty)
              val diffToChangeLog = new DiffToChangeLog(filteredDiffResult, diffOutputControl)
              diffToChangeLog.print(printStream)
              val changeSetsScoped = XML.loadString(outputStream.toString) \ "changeSet" \ "_"
              val changeSets = changeSetsScoped map stripNodeScope
              fail(changeSets.mkString(
                s"The following changes are in $schemaManager but not in $otherSchemaManager:\n  ",
                "\n  ",
                "\nEnsure that the columns/fields exist, with the same lengths in " +
                  s"$schemaManager and $otherSchemaManager and synchronize the two."))
            }
          }
        }
      }

      it should "match expected generated names" in {
        var schemaMetadata: SchemaMetadata = null

        for {
          slickDatabase <- inMemoryDatabase(databaseType, schemaManager).autoClosed
        } {
          schemaMetadata = getSchemaMetadata(slickDatabase)
        }

        var misnamed = Seq.empty[String]

        schemaMetadata.primaryKeyMetadata foreach { primaryKey =>
          val actual = primaryKey.pkName.get
          val expected = s"PK_${primaryKey.table.name}"
          if (actual != expected) {
            misnamed :+=
              s"""|  PrimaryKey: $actual
                  |  Should be:  $expected
                  |""".stripMargin
          }
        }

        schemaMetadata.foreignKeyMetadata foreach { foreignKey =>
          val actual = foreignKey.fkName.get
          val expected = s"FK_${foreignKey.fkTable.name}_${foreignKey.fkColumn}"
          if (actual != expected) {
            misnamed :+=
              s"""|  ForeignKey: $actual
                  |  Should be:  $expected
                  |""".stripMargin
          }
        }

        schemaMetadata.indexMetadata.filterNot(isForeignKeyIndex).groupBy(getIndexName) foreach {
          case (indexName, indexColumns) =>
            val index = indexColumns.head
            val prefix = if (index.nonUnique) "IX" else "UC"
            val tableName = index.table.name
            val sortedColumns = indexColumns.sortBy(_.ordinalPosition)
            val abbrColumns = sortedColumns.map(indexColumn => snakeAbbreviate(indexColumn.column.get))

            val actual = indexName
            val expected = abbrColumns.mkString(s"${prefix}_${tableName}_", "_", "")

            if (actual != expected) {
              misnamed :+=
                s"""|  Index:     $actual
                    |  Should be: $expected
                    |""".stripMargin
            }
        }

        var missing = Seq.empty[String]

        schemaMetadata.columns foreach { column =>
          if (!schemaMetadata.existsTableItem(column)) {
            missing :+= s"  ${tableClassName(column.tableName)}.${column.itemName}"
          }
        }

        schemaMetadata.slickItems foreach { databaseItem =>
          if (!schemaMetadata.existsSlickMapping(databaseItem)) {
            missing :+= s"  ${slickClassName(databaseItem.tableName)}.${databaseItem.itemName}"
          }
        }

        if (missing.nonEmpty || misnamed.nonEmpty) {
          var failMessage = ""

          if (misnamed.nonEmpty) {
            failMessage += misnamed.mkString(s"The following items are misnamed in $schemaManager:\n", "\n", "\n")
          }

          if (missing.nonEmpty) {
            failMessage += missing.mkString(
              s"Based on the schema in $schemaManager, please ensure that the following tables/columns exist:\n",
              "\n", "\n")
          }

          fail(failMessage)
        }
      }
    }
  }
}

object SchemaManagerSpec {
  // strip the namespace from elems and their children
  private def stripNodeScope(node: Node): Node = {
    node match {
      case elem: Elem => elem.copy(scope = TopScope, child = elem.child map stripNodeScope)
      case other => other
    }
  }

  private def compare[ReferenceProfile <: JdbcProfile, ComparisonProfile <: JdbcProfile, T]
  (referenceProfile: ReferenceProfile,
   referenceDatabase: ReferenceProfile#Backend#Database,
   comparisonProfile: ComparisonProfile,
   comparisonDatabase: ComparisonProfile#Backend#Database)(block: DiffResult => T): T = {
    DatabaseTestKit.withConnections(referenceProfile, referenceDatabase, comparisonProfile, comparisonDatabase) {
      LiquibaseUtils.compare(_, _)(block)
    }
  }


  private val SnakeRegex = "_([a-z])".r

  private def snakeToCamel(value: String): String = {
    SnakeRegex.replaceAllIn(value.toLowerCase, _.group(1).toUpperCase)
  }

  private def snakeAbbreviate(value: String): String = {
    SnakeRegex.findAllMatchIn("_" + value.toLowerCase).map(_.group(1)).mkString("").toUpperCase
  }

  private def tableClassName(tableName: String) = s"cromwell.database.sql.tables.$tableName"

  private def slickClassName(tableName: String) =
    s"cromwell.database.slick.tables.${tableName}Component$$${tableName.replace("Entry", "Entries")}"

  private def getIndexName(index: MIndexInfo) = index.indexName.get.replaceAll("(^SYS_IDX_|_\\d+$)", "")

  private def isForeignKeyIndex(index: MIndexInfo) = getIndexName(index).startsWith("FK_")

  case class TableClass(tableName: String) {
    private def getClass(name: String): Try[Class[_]] = Try(Class.forName(name))

    private lazy val tableColumns = getClass(tableClassName(tableName)).map(_.getDeclaredFields).getOrElse(Array.empty)
    private lazy val slickMapping = getClass(slickClassName(tableName)).map(_.getDeclaredMethods).getOrElse(Array.empty)

    def existsTableField(name: String): Boolean = tableColumns.exists(_.getName == name)

    def existsSlickMapping(name: String): Boolean = slickMapping.exists(_.getName == name)
  }

  case class DatabaseItem(tableName: String, itemName: String)

  case class SchemaMetadata(tableMetadata: Seq[MTable], columnMetadata: Seq[MColumn], indexMetadata: Seq[MIndexInfo],
                            primaryKeyMetadata: Seq[MPrimaryKey], foreignKeyMetadata: Seq[MForeignKey]) {
    lazy val tables: Seq[TableClass] = tableMetadata.map({ table =>
      val tableName = snakeToCamel(table.name.name).capitalize
      TableClass(tableName)
    }).distinct

    lazy val columns: Seq[DatabaseItem] = columnMetadata.map({ column =>
      val tableName = snakeToCamel(column.table.name).capitalize
      val columnName = snakeToCamel(column.name)
      DatabaseItem(tableName, columnName)
    }).distinct

    lazy val indexes: Seq[DatabaseItem] = indexMetadata.map({ index =>
      val tableName = snakeToCamel(index.table.name).capitalize
      val indexName = snakeToCamel(getIndexName(index))
      DatabaseItem(tableName, indexName)
    }).distinct

    lazy val foreignKeys: Seq[DatabaseItem] = foreignKeyMetadata.map({ foreignKey =>
      val tableName = snakeToCamel(foreignKey.fkTable.name).capitalize
      val indexName = snakeToCamel(foreignKey.fkName.get)
      DatabaseItem(tableName, indexName)
    }).distinct

    lazy val slickItems: Seq[DatabaseItem] = columns ++ indexes ++ foreignKeys

    def existsTableItem(tableItem: DatabaseItem): Boolean = {
      tables.find(_.tableName == tableItem.tableName).exists(_.existsTableField(tableItem.itemName))
    }

    def existsSlickMapping(tableItem: DatabaseItem): Boolean = {
      tables.find(_.tableName == tableItem.tableName).exists(_.existsSlickMapping(tableItem.itemName))
    }
  }

}
