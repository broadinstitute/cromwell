package cromwell.services.database

import cromwell.database.migration.liquibase.LiquibaseUtils
import liquibase.database.ObjectQuotingStrategy
import org.scalatest.{FlatSpec, Matchers}
import cats.data.NonEmptyList

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext

class LiquibaseChangeSetSpec extends FlatSpec with Matchers {

  implicit val executionContext = ExecutionContext.global

  behavior of "Liquibase Change Sets"

  CromwellDatabaseType.All foreach { databaseType =>
    val changeSets = LiquibaseUtils.getChangeSets(databaseType.liquibaseSettings)

    changeSets foreach { changeSet =>
      val dbmsSet = Option(changeSet.getDbmsSet).toSet[java.util.Set[String]].flatMap(_.asScala)
      val dbmsNelOption = NonEmptyList.fromList(dbmsSet.toList.sorted)

      val description = s"${databaseType.name} change set" +
        s" ${changeSet.getId} in ${changeSet.getChangeLog} by ${changeSet.getAuthor}" +
        s" for ${dbmsNelOption.map(_.toList.mkString("(", ",", ")")).getOrElse("[dbms not set]")}"

      it should s"check dbms in $description" in {
        withClue("the dbms attribute must explicitly list supported databases:") {
          changeSet.getDbmsSet shouldNot be(null)
          changeSet.getDbmsSet shouldNot be(empty)
          changeSet.getDbmsSet shouldNot contain atLeastOneOf("none", "all")
        }
        changeSet.getDbmsSet.asScala foreach { dbms =>
          withClue("do not use dbms excludes:") {
            dbms shouldNot include("!")
          }
        }
      }

      if (dbmsSet.contains("postgresql")) {
        it should s"check object quoting in $description" in {
          withClue("""databaseChangeLog objectQuotingStrategy="QUOTE_ALL_OBJECTS" must be set for PostgreSQL:""") {
            changeSet.getObjectQuotingStrategy should be (ObjectQuotingStrategy.QUOTE_ALL_OBJECTS)
          }
        }
      }
    }
  }
}
