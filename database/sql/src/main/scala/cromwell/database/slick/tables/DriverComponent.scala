package cromwell.database.slick.tables

import java.sql.{Blob, Clob}

import javax.sql.rowset.serial.{SerialBlob, SerialClob}
import org.apache.commons.io.IOUtils
import slick.jdbc.{JdbcProfile, MySQLProfile, PostgresProfile}

trait DriverComponent {
  val driver: JdbcProfile

  import driver.api._

  /** Ensure clobs are retrieved inside the transaction, not after */
  implicit val serialClobColumnType = MappedColumnType.base[SerialClob, Clob](
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
  )

  /** Ensure clobs are retrieved inside the transaction, not after */
  implicit val serialBlobColumnType = MappedColumnType.base[SerialBlob, Blob](
    identity,
    {
      case serialBlob: SerialBlob => serialBlob
      case blob => new SerialBlob(blob)
    }
  )

  private val shouldQuote = this.driver match {
    // https://stackoverflow.com/questions/43111996/why-postgresql-does-not-like-uppercase-table-names#answer-43112096
    case PostgresProfile => true
    case _ => false
  }

  /** Adds quotes around the string if required by the DBMS. */
  def quoted(string: String) = if (shouldQuote) s""""$string"""" else string

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

}
