package cromwell.services.database

import cromwell.core.Tags._
import cromwell.database.slick.SlickDatabase
import cromwell.services.database.LiquibaseComparisonSpec._
import cromwell.services.database.LiquibaseOrdering._
import liquibase.snapshot.DatabaseSnapshot
import liquibase.statement.DatabaseFunction
import liquibase.structure.DatabaseObject
import liquibase.structure.core._
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}
import slick.jdbc.GetResult

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.reflect._

/**
  * Compares all of the various liquibase schemas against an in-memory HSQLDB-Slick schema.
  */
class LiquibaseComparisonSpec extends FlatSpec with Matchers with ScalaFutures {

  implicit val executionContext = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = scaled(5.seconds), interval = scaled(100.millis))

  CromwellDatabaseType.All foreach { databaseType =>
    lazy val expectedSnapshot = DatabaseTestKit.inMemorySnapshot(databaseType, SlickSchemaManager)
    lazy val expectedColumns = get[Column](expectedSnapshot)
    lazy val expectedPrimaryKeys = get[PrimaryKey](expectedSnapshot)
    lazy val expectedForeignKeys = get[ForeignKey](expectedSnapshot)
    lazy val expectedUniqueConstraints = get[UniqueConstraint](expectedSnapshot)
    lazy val expectedIndexes = get[Index](expectedSnapshot) filterNot DatabaseTestKit.isGenerated

    DatabaseSystem.All foreach { databaseSystem =>

      behavior of s"Liquibase Comparison for ${databaseType.name} ${databaseSystem.shortName}"

      lazy val liquibasedDatabase = DatabaseTestKit.initializedDatabaseFromSystem(databaseType, databaseSystem)

      lazy val actualSnapshot = DatabaseTestKit.liquibaseSnapshot(liquibasedDatabase)
      lazy val actualColumns = get[Column](actualSnapshot)
      lazy val actualPrimaryKeys = get[PrimaryKey](actualSnapshot)
      lazy val actualForeignKeys = get[ForeignKey](actualSnapshot)
      lazy val actualUniqueConstraints = get[UniqueConstraint](actualSnapshot)
      lazy val actualIndexes = get[Index](actualSnapshot)

      lazy val columnMapping = getColumnMapping(databaseSystem)

      expectedColumns foreach { expectedColumn =>
        val description = s"column ${expectedColumn.getRelation.getName}.${expectedColumn.getName}"

        it should s"match the Slick schema for $description" taggedAs DbmsTest in {
          val actualColumnOption = actualColumns find { actualColumn =>
            ColumnDescription.from(actualColumn) == ColumnDescription.from(expectedColumn)
          }
          val actualColumn = actualColumnOption getOrElse fail(s"Did not find $description")

          withClue(s"for type " +
            s"${actualColumn.getType.getTypeName}(default = ${actualColumn.getDefaultValue}) vs. " +
            s"${expectedColumn.getType.getTypeName}(default = ${expectedColumn.getDefaultValue}):") {

            val actualColumnType = ColumnType.from(actualColumn)

            actualColumn.isAutoIncrement should be(expectedColumn.isAutoIncrement)
            if (expectedColumn.isAutoIncrement) {
              // Auto increment columns may have different types, such as SERIAL/BIGSERIAL
              // https://www.postgresql.org/docs/11/datatype-numeric.html#DATATYPE-SERIAL
              val actualColumnDefault = ColumnDefault(actualColumnType, actualColumn.getDefaultValue)
              val autoIncrementDefault = getAutoIncrementDefault(databaseSystem, columnMapping, expectedColumn)
              actualColumnDefault should be(autoIncrementDefault)
            } else {

              // Check the column type
              val mappedColumnType = getColumnType(expectedColumn, columnMapping)
              actualColumnType should be(mappedColumnType)

              if (isSlickDefaultNull(expectedColumn)) {
                // The column has a default value of "null"
                // There are a number of ways off expressing null/None/NULL when dealing with CLOB/BLOB/BIT etc.
                val expectedOptions = Set(
                  None,
                  Option(DefaultNullBoolean),
                  Option(DefaultNullString),
                  Option(DefaultNullFunction),
                )
                List(Option(actualColumn.getDefaultValue)) should contain atLeastOneElementOf expectedOptions
              } else {

                // Check the default value
                val mappedDefaultValue = getColumnDefault(expectedColumn, columnMapping)
                actualColumn.getDefaultValue should be(mappedDefaultValue)
              }
            }

            // Check for column nullability
            val nullTodos = getNullTodos(databaseSystem, databaseType)
            if (nullTodos.contains(ColumnDescription.from(actualColumn))) {
              // Oops. This column is nullable. TODO: make a changelog to fix, and then remove it from the list
              assert(actualColumn.isNullable, "Column is explicitly listed as a null field:")
            } else {
              actualColumn.isNullable should be(expectedColumn.isNullable)
            }

            // Some types aren't available via liquibase, so go into the database to get the types
            columnTypeValidationOption(expectedColumn, databaseSystem) foreach { expectedColumnType =>
              val dbio = columnTypeDbio(expectedColumn, databaseSystem, liquibasedDatabase)
              val future = liquibasedDatabase.database.run(dbio)
              future.futureValue should be(expectedColumnType)
            }

            // Verify that sequence widths are the same as columns
            sequenceTypeValidationOption(expectedColumn, databaseSystem) foreach { expectedSequenceType =>
              val dbio = sequenceTypeDbio(expectedColumn, databaseSystem, liquibasedDatabase)
              val future = liquibasedDatabase.database.run(dbio)
              val res = future.futureValue
              res should be(expectedSequenceType)
            }

          }
        }
      }

      expectedPrimaryKeys foreach { expectedPrimaryKey =>
        val description = s"primary key on ${expectedPrimaryKey.getTable.getName}"

        it should s"match the Slick schema for $description" taggedAs DbmsTest in {
          val actualPrimaryKeyOption = actualPrimaryKeys find {
            _.getTable.getName == expectedPrimaryKey.getTable.getName
          }
          val actualPrimaryKey = actualPrimaryKeyOption getOrElse fail(s"Did not find $description")

          actualPrimaryKey.getColumns.asScala.map(ColumnDescription.from) should
            contain theSameElementsAs expectedPrimaryKey.getColumns.asScala.map(ColumnDescription.from)
        }
      }

      expectedForeignKeys foreach { expectedForeignKey =>
        val description = s"foreign key ${expectedForeignKey.getName}"

        it should s"match the Slick schema for $description" taggedAs DbmsTest in {
          val actualForeignKeyOption = actualForeignKeys find {
            _.getName == expectedForeignKey.getName
          }
          val actualForeignKey = actualForeignKeyOption getOrElse fail(s"Did not find $description")

          actualForeignKey.getPrimaryKeyColumns.asScala.map(ColumnDescription.from) should
            contain theSameElementsAs expectedForeignKey.getPrimaryKeyColumns.asScala.map(ColumnDescription.from)
          actualForeignKey.getForeignKeyColumns.asScala.map(ColumnDescription.from) should
            contain theSameElementsAs expectedForeignKey.getForeignKeyColumns.asScala.map(ColumnDescription.from)

          expectedForeignKey.getUpdateRule match {
            case ForeignKeyConstraintType.importedKeyRestrict | ForeignKeyConstraintType.importedKeyNoAction =>
              actualForeignKey.getUpdateRule should
                (be(ForeignKeyConstraintType.importedKeyRestrict) or be(ForeignKeyConstraintType.importedKeyNoAction))
            case other => actualForeignKey.getUpdateRule should be(other)
          }
          expectedForeignKey.getDeleteRule match {
            case ForeignKeyConstraintType.importedKeyRestrict | ForeignKeyConstraintType.importedKeyNoAction =>
              actualForeignKey.getDeleteRule should
                (be(ForeignKeyConstraintType.importedKeyRestrict) or be(ForeignKeyConstraintType.importedKeyNoAction))
            case other => actualForeignKey.getDeleteRule should be(other)
          }
          actualForeignKey.isDeferrable should be(expectedForeignKey.isDeferrable)
          actualForeignKey.isInitiallyDeferred should be(expectedForeignKey.isInitiallyDeferred)
        }
      }

      expectedUniqueConstraints foreach { expectedUniqueConstraint =>
        val description = s"unique constraint ${expectedUniqueConstraint.getName}"

        it should s"match the Slick schema for $description" taggedAs DbmsTest in {
          val actualUniqueConstraintOption = actualUniqueConstraints find {
            _.getName == expectedUniqueConstraint.getName
          }
          val actualUniqueConstraint = actualUniqueConstraintOption getOrElse
            fail(s"Did not find $description")

          actualUniqueConstraint.getColumns.asScala.map(ColumnDescription.from) should
            contain theSameElementsAs expectedUniqueConstraint.getColumns.asScala.map(ColumnDescription.from)
        }
      }

      expectedIndexes foreach { expectedIndex =>
        val description = s"index ${expectedIndex.getName}"

        it should s"match the Slick schema for $description" taggedAs DbmsTest in {
          val actualIndexOption = actualIndexes find {
            _.getName == expectedIndex.getName
          }
          val actualIndex = actualIndexOption getOrElse fail(s"Did not find $description")

          actualIndex.isUnique should be(expectedIndex.isUnique)
          actualIndex.getColumns.asScala.map(ColumnDescription.from) should
            contain theSameElementsAs expectedIndex.getColumns.asScala.map(ColumnDescription.from)
        }
      }

      it should "close the database" taggedAs DbmsTest in {
        liquibasedDatabase.close()
      }
    }
  }
}

object LiquibaseComparisonSpec {
  private def get[T <: DatabaseObject : ClassTag : Ordering](databaseSnapshot: DatabaseSnapshot): Seq[T] = {
    val databaseObjectClass = classTag[T].runtimeClass.asInstanceOf[Class[T]]
    databaseSnapshot.get(databaseObjectClass).asScala.toSeq.sorted
  }

  private val DefaultNullBoolean = Boolean.box(false)
  private val DefaultNullString = "NULL"
  private val DefaultNullFunction = new DatabaseFunction(DefaultNullString)

  private def isSlickDefaultNull(column: Column): Boolean = {
    Option(column.getDefaultValue).isEmpty || column.getDefaultValue == DefaultNullFunction
  }

  case class ColumnDescription(tableName: String, columnName: String)

  object ColumnDescription {
    def from(column: Column): ColumnDescription = {
      ColumnDescription(column.getRelation.getName, column.getName)
    }
  }

  case class ColumnType
  (
    typeName: String,
    sizeOption: Option[Int] = None,
  )

  object ColumnType {
    def from(column: Column): ColumnType = {
      ColumnType(
        column.getType.getTypeName.toUpperCase,
        Option(column.getType.getColumnSize).map(_.toInt),
      )
    }
  }

  case class ColumnDefault
  (
    columnType: ColumnType,
    defaultValue: AnyRef,
  )

  object ColumnDefault {
    def from(column: Column): ColumnDefault = {
      ColumnDefault(ColumnType.from(column), column.getDefaultValue)
    }
  }

  case class ColumnMapping
  (
    typeMapping: Map[ColumnType, ColumnType] = Map.empty,
    defaultMapping: Map[ColumnDefault, ColumnDefault] = Map.empty,
  )

  /** Generate the expected PostgreSQL sequence name for a column. */
  private def postgresqlSeqName(column: Column): String = {

    def pad(name: String): String = if (name.endsWith("_")) name else name + "_"

    // Postgres cuts of the length of names around this length
    val Count = 30

    def shorten(name: String, isColumn: Boolean): String = {
      pad {
        // NOTE: Table and column name truncation seems slightly different.
        // This logic was empirically derived. Feel free to modify/simplify!
        if (name.length < Count || isColumn && name.length == Count) {
          name
        } else if (isColumn && name.length == Count + 1) {
          name.dropRight(1)
        } else {
          name.take(Count - 1)
        }
      }
    }

    val tableName = shorten(column.getRelation.getName, isColumn = false)
    val columnName = shorten(column.getName, isColumn = true)
    s"$tableName${columnName}seq"
  }

  // Types as they are represented in HSQLDB that will have different representations in other DBMS.
  private val HsqldbTypeBigInt = ColumnType("BIGINT", Option(64))
  private val HsqldbTypeBlob = ColumnType("BLOB", Option(1073741824))
  private val HsqldbTypeBoolean = ColumnType("BOOLEAN", Option(0))
  private val HsqldbTypeClob = ColumnType("CLOB", Option(1073741824))
  private val HsqldbTypeInteger = ColumnType("INTEGER", Option(32))
  private val HsqldbTypeTimestamp = ColumnType("TIMESTAMP")

  // Nothing to map as the original is also HSQLDB
  private val HsqldbColumnMapping = ColumnMapping()

  // Note: BIT vs. TINYINT may be yet another tabs vs. spaces
  // https://stackoverflow.com/questions/11167793/boolean-or-tinyint-confusion/17298805
  private val MysqldbColumnMapping =
  ColumnMapping(
    typeMapping = Map(
      HsqldbTypeBigInt -> ColumnType("BIGINT", Option(19)),
      HsqldbTypeBlob -> ColumnType("LONGBLOB", Option(2147483647)),
      HsqldbTypeBoolean -> ColumnType("TINYINT", Option(3)),
      HsqldbTypeClob -> ColumnType("LONGTEXT", Option(2147483647)),
      HsqldbTypeInteger -> ColumnType("INT", Option(10)),
      HsqldbTypeTimestamp -> ColumnType("DATETIME"),
    ),
    defaultMapping = Map(
      ColumnDefault(HsqldbTypeBoolean, Boolean.box(true)) ->
        ColumnDefault(ColumnType("TINYINT", Option(3)), Int.box(1)),
    ),
  )

  // MariaDB should behave exactly the same as MySQL
  private val MariadbColumnMapping = MysqldbColumnMapping

  private val PostgresqlColumnMapping =
    ColumnMapping(
      typeMapping = Map(
        HsqldbTypeBigInt -> ColumnType("INT8", None),
        HsqldbTypeBlob -> ColumnType("OID", None),
        HsqldbTypeBoolean -> ColumnType("BOOL", None),
        HsqldbTypeClob -> ColumnType("TEXT", None),
        HsqldbTypeInteger -> ColumnType("INT4", None),
      ),
    )

  /**
    * Returns the column mapping for the DBMS.
    */
  private def getColumnMapping(databaseSystem: DatabaseSystem): ColumnMapping = {
    databaseSystem match {
      case HsqldbDatabaseSystem => HsqldbColumnMapping
      case MariadbDatabaseSystem => MariadbColumnMapping
      case MysqlDatabaseSystem => MysqldbColumnMapping
      case PostgresqlDatabaseSystem => PostgresqlColumnMapping
    }
  }

  /**
    * Returns the column type, possibly mapped via the ColumnMapping.
    */
  private def getColumnType(column: Column, columnMapping: ColumnMapping): ColumnType = {
    val columnType = ColumnType.from(column)
    columnMapping.typeMapping.getOrElse(columnType, columnType)
  }

  /**
    * Returns the default for the column, either from ColumnMapping or the column itself.
    */
  private def getColumnDefault(column: Column, columnMapping: ColumnMapping): AnyRef = {
    columnMapping.defaultMapping get ColumnDefault.from(column) map (_.defaultValue) getOrElse column.getDefaultValue
  }

  /**
    * Return the default for the auto increment column.
    */
  private def getAutoIncrementDefault(databaseSystem: DatabaseSystem,
                                      columnMapping: ColumnMapping,
                                      column: Column): ColumnDefault = {
    databaseSystem match {
      case PostgresqlDatabaseSystem =>
        val columnType = column.getType.getTypeName match {
          case "BIGINT" => ColumnType("BIGSERIAL", None)
          case "INTEGER" => ColumnType("SERIAL", None)
        }
        val columnDefault = new DatabaseFunction(s"""nextval('"${postgresqlSeqName(column)}"'::regclass)""")
        ColumnDefault(columnType, columnDefault)
      case _ => ColumnDefault(getColumnType(column, columnMapping), column.getDefaultValue)
    }
  }

  /**
    * Returns an optional extra check to ensure that datetimes can store microseconds.
    *
    * This is to double check for MySQL:
    *
    * > ... the default precision is 0.
    * > This differs from the standard SQL default of 6, for compatibility with previous MySQL versions.
    *
    * via: https://dev.mysql.com/doc/refman/8.0/en/fractional-seconds.html
    *
    * And for MariaDB:
    *
    * > The microsecond precision can be from 0-6. If not specified 0 is used.
    *
    * via: https://mariadb.com/kb/en/library/datetime/
    *
    * This check also has to be done here, as Liquibase does not return the precision for Mysql datetime fields.
    */
  private def columnTypeValidationOption(column: Column, databaseSystem: DatabaseSystem): Option[String] = {
    databaseSystem match {
      case MysqlDatabaseSystem | MariadbDatabaseSystem if column.getType.getTypeName == "TIMESTAMP" =>
        Option("datetime(6)")
      case _ => None
    }
  }

  private def columnTypeDbio(column: Column,
                             databaseSystem: DatabaseSystem,
                             database: SlickDatabase): database.dataAccess.driver.api.DBIO[String] = {
    import database.dataAccess.driver.api._
    databaseSystem match {
      case MysqlDatabaseSystem | MariadbDatabaseSystem if column.getType.getTypeName == "TIMESTAMP" =>
        val getType = GetResult(_.rs.getString("Type"))

        //noinspection SqlDialectInspection
        sql"""SHOW COLUMNS
              FROM #${column.getRelation.getName}
              WHERE FIELD = '#${column.getName}'
           """.as[String](getType).head
      case _ => DBIO.failed(unsupportedColumnTypeException(column, databaseSystem))
    }
  }

  /**
    * Returns an optional extra check to ensure that sequences have the same types as their auto increment columns.
    *
    * This is because PostgreSQL requires two statements to modify SERIAL columns to BIGSERIAL, one to widen the column,
    * and another to widen the sequence.
    *
    * https://stackoverflow.com/questions/52195303/postgresql-primary-key-id-datatype-from-serial-to-bigserial#answer-52195920
    * https://www.postgresql.org/docs/11/datatype-numeric.html#DATATYPE-SERIAL
    */
  private def sequenceTypeValidationOption(column: Column, databaseSystem: DatabaseSystem): Option[String] = {
    databaseSystem match {
      case PostgresqlDatabaseSystem if column.isAutoIncrement => Option(column.getType.getTypeName.toLowerCase)
      case _ => None
    }
  }

  private def sequenceTypeDbio(column: Column,
                               databaseSystem: DatabaseSystem,
                               database: SlickDatabase): database.dataAccess.driver.api.DBIO[String] = {
    import database.dataAccess.driver.api._
    databaseSystem match {
      case PostgresqlDatabaseSystem if column.isAutoIncrement =>

        //noinspection SqlDialectInspection
        sql"""select data_type
              from INFORMATION_SCHEMA.sequences
              where sequence_name = '#${postgresqlSeqName(column)}'
           """.as[String].head
      case _ => DBIO.failed(unsupportedColumnTypeException(column, databaseSystem))
    }
  }

  private def unsupportedColumnTypeException(column: Column,
                                             databaseSystem: DatabaseSystem): UnsupportedOperationException = {
    new UnsupportedOperationException(
      s"${databaseSystem.shortName} ${column.getRelation.getName}.${column.getName}: ${column.getType.getTypeName}"
    )
  }

  /**
    * Returns columns that are nullable, but shouldn't be.
    *
    * These nullables are not really hurting anything, but we should not add any more columns to this list.
    *
    * TODO: make a changelog to fix, and then remove list of mistakes.
    */
  private def getNullTodos(databaseSystem: DatabaseSystem,
                           databaseType: CromwellDatabaseType[_ <: SlickDatabase]): Seq[ColumnDescription] = {
    (databaseSystem, databaseType) match {
      case (MysqlDatabaseSystem, EngineDatabaseType) =>
        List(
          ColumnDescription("CALL_CACHING_DETRITUS_ENTRY", "CALL_CACHING_ENTRY_ID"),
          ColumnDescription("CALL_CACHING_DETRITUS_ENTRY", "DETRITUS_KEY"),
          ColumnDescription("CALL_CACHING_ENTRY", "ALLOW_RESULT_REUSE"),
          ColumnDescription("CALL_CACHING_ENTRY", "CALL_FULLY_QUALIFIED_NAME"),
          ColumnDescription("CALL_CACHING_ENTRY", "JOB_ATTEMPT"),
          ColumnDescription("CALL_CACHING_ENTRY", "JOB_INDEX"),
          ColumnDescription("CALL_CACHING_ENTRY", "WORKFLOW_EXECUTION_UUID"),
          ColumnDescription("CALL_CACHING_HASH_ENTRY", "CALL_CACHING_ENTRY_ID"),
          ColumnDescription("CALL_CACHING_SIMPLETON_ENTRY", "CALL_CACHING_ENTRY_ID"),
          ColumnDescription("DOCKER_HASH_STORE_ENTRY", "DOCKER_SIZE"),
          ColumnDescription("JOB_KEY_VALUE_ENTRY", "CALL_FULLY_QUALIFIED_NAME"),
          ColumnDescription("JOB_KEY_VALUE_ENTRY", "JOB_ATTEMPT"),
          ColumnDescription("JOB_KEY_VALUE_ENTRY", "JOB_INDEX"),
          ColumnDescription("DOCKER_HASH_STORE_ENTRY", "DOCKER_SIZE"),
          ColumnDescription("JOB_STORE_ENTRY", "CALL_FULLY_QUALIFIED_NAME"),
          ColumnDescription("JOB_STORE_ENTRY", "JOB_ATTEMPT"),
          ColumnDescription("JOB_STORE_ENTRY", "JOB_INDEX"),
          ColumnDescription("JOB_STORE_ENTRY", "RETRYABLE_FAILURE"),
          ColumnDescription("JOB_STORE_ENTRY", "WORKFLOW_EXECUTION_UUID"),
          ColumnDescription("JOB_STORE_SIMPLETON_ENTRY", "JOB_STORE_ENTRY_ID"),
          ColumnDescription("WORKFLOW_STORE_ENTRY", "IMPORTS_ZIP"),
          ColumnDescription("WORKFLOW_STORE_ENTRY", "WORKFLOW_EXECUTION_UUID"),
          ColumnDescription("WORKFLOW_STORE_ENTRY", "WORKFLOW_STATE"),
        )
      case (MysqlDatabaseSystem, MetadataDatabaseType) =>
        List(
          ColumnDescription("CUSTOM_LABEL_ENTRY", "CUSTOM_LABEL_KEY"),
          ColumnDescription("CUSTOM_LABEL_ENTRY", "CUSTOM_LABEL_VALUE"),
          ColumnDescription("SUMMARY_STATUS_ENTRY", "SUMMARY_NAME"),
          ColumnDescription("SUMMARY_STATUS_ENTRY", "SUMMARY_POSITION"),
        )
      case _ => Nil
    }
  }
}
