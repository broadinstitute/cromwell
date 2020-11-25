package cromwell.database.slick.tables

import java.sql._

import cromwell.database.sql.SqlConverters._
import javax.sql.rowset.serial.{SerialBlob, SerialClob}
import org.apache.commons.io.IOUtils
import slick.jdbc.{JdbcProfile, MySQLProfile, PostgresProfile, SQLiteProfile}

import scala.language.higherKinds
import scala.reflect.ClassTag

trait DriverComponent {
  driverComponent =>
  val driver: JdbcProfile

  import driver.api._

  /**
    * Use String instead of CLOB on SQLite.
    * https://github.com/xerial/sqlite-jdbc/blob/3.32.3.2/src/main/java/org/sqlite/jdbc3/JDBC3PreparedStatement.java#L513-L514
    */
  implicit val clobOrStringColumnType: driver.DriverJdbcType[Clob] =
    new driver.columnTypes.ClobJdbcType {
      override def setValue(v: Clob, p: PreparedStatement, idx: Int): Unit = {
        driverComponent.driver match {
          case SQLiteProfile => p.setString(idx, v.toRawString)
          case _ => super.setValue(v, p, idx)
        }
      }

      override def getValue(r: ResultSet, idx: Int): Clob = {
        driverComponent.driver match {
          case SQLiteProfile => r.getString(idx).toClobOption.orNull
          case _ => super.getValue(r, idx)
        }
      }

      override def updateValue(v: Clob, r: ResultSet, idx: Int): Unit = {
        driverComponent.driver match {
          case SQLiteProfile => r.updateString(idx, v.toRawString)
          case _ => super.updateValue(v, r, idx)
        }
      }
    }

  /**
    * Use Array[Byte] instead of BLOB on SQLite.
    * https://github.com/xerial/sqlite-jdbc/blob/3.32.3.2/src/main/java/org/sqlite/jdbc3/JDBC3PreparedStatement.java#L511-L512
    */
  implicit val blobOrBytesColumnType: driver.DriverJdbcType[Blob] =
    new driver.columnTypes.BlobJdbcType {
      override def setValue(v: Blob, p: PreparedStatement, idx: Int): Unit = {
        driverComponent.driver match {
          case SQLiteProfile => p.setBytes(idx, v.toBytes)
          case _ => super.setValue(v, p, idx)
        }
      }

      override def getValue(r: ResultSet, idx: Int): Blob = {
        driverComponent.driver match {
          case SQLiteProfile => r.getBytes(idx).toBlobOption.orNull
          case _ => super.getValue(r, idx)
        }
      }

      override def updateValue(v: Blob, r: ResultSet, idx: Int): Unit = {
        driverComponent.driver match {
          case SQLiteProfile => r.updateBytes(idx, v.toBytes)
          case _ => super.updateValue(v, r, idx)
        }
      }
    }

  /** Ensure clobs are retrieved inside the transaction, not after */
  implicit val serialClobColumnType: BaseColumnType[SerialClob] = MappedColumnType.base[SerialClob, Clob](
    identity,
    {
      case serialClob: SerialClob => serialClob
      case clob =>
        /*
        PostgreSQL's JDBC driver has issues with non-ascii characters.
        https://stackoverflow.com/questions/5043992/postgres-utf-8-clobs-with-jdbc

        It returns bad values for length() and getAsciiStream(), and causes an extra null bytes to be added at the end
        of the resultant SerialClob.

        Example via copy_workflow_outputs/unscattered.wdl:

          "... Enfin un peu de francais pour contrer ce raz-de-marée anglais ! ..."

        The 'é' in results in an extra null byte at the end of getAsciiStream().
         */
        val string = IOUtils.toString(clob.getCharacterStream)
        new SerialClob(string.toCharArray)
    }
  )(ClassTag(classOf[SerialClob]), clobOrStringColumnType)

  /** Ensure blobs are retrieved inside the transaction, not after */
  implicit val serialBlobColumnType: BaseColumnType[SerialBlob] = MappedColumnType.base[SerialBlob, Blob](
    identity,
    {
      case serialBlob: SerialBlob => serialBlob
      case blob => new SerialBlob(blob)
    }
  )(ClassTag(classOf[SerialBlob]), blobOrBytesColumnType)

  private val shouldQuote = this.driver match {
    // https://stackoverflow.com/questions/43111996/why-postgresql-does-not-like-uppercase-table-names#answer-43112096
    case PostgresProfile => true
    case _ => false
  }

  /** Adds quotes around the string if required by the DBMS. */
  def quoted(string: String): String = if (shouldQuote) s""""$string"""" else string

  val clobToString: Rep[SerialClob] => Rep[String] = {
    this.driver match {
        /*
        Workaround https://jira.mariadb.org/browse/CONJ-717
        Bypass Slick `asColumnOf[String]` calling the JDBC `{fn convert(Column, VARCHAR)}`.
        Instead directly call `concat(Column)` supported by both the MariaDB driver and the MySQL driver.
         */
      case MySQLProfile => SimpleFunction.unary("concat")
      case _ => _.asColumnOf[String]
    }
  }

  /**
    * Disable "SELECT ... FOR UPDATE" on SQLite as it doesn't support locking specifically for updates:
    * https://stackoverflow.com/questions/46945790/why-could-select-for-update-simply-be-ignored-with-databases-that-dont-support
    */
  def maybeSelectForUpdate[E, U, C[_]](query: Query[E, U, C]): Query[E, U, C] = {
    this.driver match {
      case SQLiteProfile => query
      case _ => query.forUpdate
    }
  }
}
