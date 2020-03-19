package cromwell

import java.time.OffsetDateTime
import java.time.format.DateTimeFormatter

import org.apache.commons.math3.fitting.WeightedObservedPoints
import org.scalatest.FlatSpec
import Helper._

class DatabaseFullnessSpec extends FlatSpec {
  behavior of "Production Database Fullness"
  // This assumes growth of the database can be described by a function of the form y = a + b*e^(c*x),
  // where y is the database size and x is time.
  // Apache Commons Math does not support exponential curve fitting directly but luckily there are maths
  // (tl;dr take logs and do linear regression):
  // https://stackoverflow.com/questions/41783768/determine-coefficients-of-exponential-function-using-apache-common-math

  it should "be able to warn about impending database explosion" ignore {
    val googleProject = System.getenv("GOOGLE_PROJECT")
    val cloudSqlInstance = System.getenv("CLOUD_SQL_INSTANCE")
    val maximumDatabaseSizeTiB = System.getenv("MAXIMUM_DATABASE_SIZE_TIB").toLong
    val databaseWatchIntervalMonths = System.getenv("DATABASE_WATCH_INTERVAL_MONTHS").toLong
    val now = OffsetDateTime.now()

    val curlCommand = buildCurlCommand(now, googleProject, cloudSqlInstance)
    // Drill into the time series points and extract the timestamp and size, format as TSV rows.
    val jqCommand = "jq -r '.timeSeries[0].points | .[] | [.interval.startTime, .value.doubleValue] | @tsv'"

    val command = List(
      "/bin/sh",
      "-c",
      s"$curlCommand | $jqCommand"
    )

    val stdout = new NewlineAddingAppender()
    val stderr = new NewlineAddingAppender()

    import sys.process._

    val rc = command ! ProcessLogger(stdout.append, stderr.append)
    if (rc != 0) {
      val message =
        s"""
           |rc is $rc
           |curlCommand is $curlCommand
           |stdout is $stdout
           |stderr is $stderr
           |""".stripMargin.trim
      println(message)
      System.exit(1)
    }

    val timestampAndByteses: Array[TimestampAndBytes] = stdout.stringBuilder.toString().split("\n") map { line =>
      line.split("\t").toList match {
        case timestamp :: bytes :: Nil =>
          TimestampAndBytes(timestamp, bytes.toDouble)
        case x =>
          throw new RuntimeException(s"GAHHHH omg an $x")
      }
    }
    import org.apache.commons.math3.fitting._
    val points: WeightedObservedPoints = new WeightedObservedPoints();
    timestampAndByteses foreach { _.addTo(points, now) }

    val fitter = PolynomialCurveFitter.create(1)
    val Array(logB, c) = fitter.fit(points.toList)
    val b = Math.exp(logB)
    val a = b - timestampAndByteses(0).bytes.inTiB
    println(s"solving for a, b, c in <size in TiB> = a + b * exp(c * <months from now>): a = $a, b = $b, c = $c")

    val sizeWatchIntervalMonthsFromNow = a + b * Math.exp(c * databaseWatchIntervalMonths)
    // y = <db size in TiB>, x = <months from now>
    // To find out when the database will explode, set y to the maximum size and solve for x.
    // y = a + b * e^(c*x)          // original
    // y - a = b * e^(c*x)          // subtract a from both sides
    // ln(y - a) = ln(b) + c*x      // take logarithms of both sides
    // ln(y - a) - ln(b) = c*x      // subtract ln(b) from both sides
    // (ln(y - a) - ln(b))/c = x    // divide both sides by c
    // x = (ln(y - a) - ln(b))/c    // swap sides
    val whenExplodeMonths = (Math.log(maximumDatabaseSizeTiB - a) - Math.log(b)) / c

    val explodeMessage = f"""
      |At recent growth rates, in $databaseWatchIntervalMonths months the database size should be about $sizeWatchIntervalMonthsFromNow%.2f TiB.
      |The database '$cloudSqlInstance' is predicted to exceed its maximum size of $maximumDatabaseSizeTiB TiB in $whenExplodeMonths%.2f months. Please plan accordingly.
    """.stripMargin
    println(explodeMessage)

    if (sizeWatchIntervalMonthsFromNow > maximumDatabaseSizeTiB || whenExplodeMonths < databaseWatchIntervalMonths) {
      System.exit(1)
    }
  }

  def buildCurlCommand(now: OffsetDateTime, googleProject: String, cloudSqlInstance: String): String = {
    // The time series endpoint is used by the CloudSQL database graphs and retains at least 30 days worth of data but
    // not much more than that. Asking for 60 days worth of data seems to get 40-50 days.
    val maxDaysAgo = 60L
    val formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd")
    val startTime = OffsetDateTime.now().minusDays(maxDaysAgo).format(formatter);
    val endTime = now.format(formatter)
    s"""curl --silent --location --header "Authorization: Bearer $$(gcloud auth print-access-token)" 'https://monitoring.clients6.google.com/v3/projects/${googleProject}/timeSeries?filter=metric.type%3D%22cloudsql.googleapis.com%2Fdatabase%2Fdisk%2Fbytes_used%22%20AND%20(resource.label.database_id%3D%22${googleProject}%3A${cloudSqlInstance}%22)&interval.startTime=${startTime}T18%3A00%3A00Z&interval.endTime=${endTime}T18%3A00%3A00Z&aggregation.alignmentPeriod=%2B10800s&aggregation.perSeriesAligner=ALIGN_MEAN'"""
  }
}

object Helper {
  implicit class EnhancedDouble(val double: Double) extends AnyVal {
    def inTiB: Double = double / (1L << 40)
  }

  case class TimestampAndBytes(timestamp: String, bytes: Double) {
    lazy val epochSeconds = OffsetDateTime.parse(timestamp).toEpochSecond

    private def yearsSince(_then: OffsetDateTime): Double = {
      // 365 = days / year
      //  24 = hours / day
      //  60 = minutes / hour
      //  60 = seconds / minute
      // multiplying through gives seconds / year
      // dividing seconds by this gives years
      (epochSeconds - _then.toEpochSecond) / (365 * 24 * 60 * 60.0)
    }

    private def monthsSince(_then: OffsetDateTime): Double = yearsSince(_then) * 12

    def addTo(points: WeightedObservedPoints, now: OffsetDateTime): Unit = {
      // The underlying growth function is assumed to be of the form y = a + b*e^(c*x). For curve fitting, initially
      // set a to zero, later on set time to zero to solve for a. After taking logarithms of both sides:
      // ln(y) = ln(b) + c*x
      // So take logarithms of the y values (database size).
      // Months and TiB seem the most natural units to describe the growth of this db.
      points.add(monthsSince(now), Math.log(bytes.inTiB))
    }
  }

  class NewlineAddingAppender() {
    val stringBuilder = new StringBuilder()
    def append(s: String) = {
      stringBuilder append s; stringBuilder append "\n"; ()
    }
  }
}
