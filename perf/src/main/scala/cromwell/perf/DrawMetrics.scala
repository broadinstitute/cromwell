package cromwell.perf

import java.nio.file.Paths
import java.util.concurrent.Executors

import cats.effect.{ExitCode, IO, IOApp}
import cats.implicits._
import com.typesafe.scalalogging.StrictLogging
import fs2.{io, text}
import scalax.chart.XYChart

import scala.concurrent.ExecutionContext

/**
  * Utility class to draw results from statsd logs coming from perf tests.
  * Arguments: metric_to_filter output_path.png statsd.log...
  * The filter is used as a "does the metric path contains that string"
  * The third to last arguments should be paths to statsd log files. They will all be plotted on the same chart.
  * e.g arguments: metadata.GET.200.p50 metadata-response-median.png statsd-1.log statsd-2.log
  */
object DrawMetrics extends IOApp with scalax.chart.module.Charting with StrictLogging {
  private val MB = 1024 * 1024
  
  def run(args: List[String]): IO[ExitCode] = {
    implicit val ec = ExecutionContext.fromExecutorService(Executors.newFixedThreadPool(5))
    implicit val contextShift = IO.contextShift(ec)

    val filter = args.head
    val output = args.drop(1).head
    val sources = args.drop(2)

    sources
      .parTraverse[IO, IO.Par, DrawMetrics.XYSeries](makeXYSeries(filter))
      .map(XYLineChart(_))
      .map(configureChart(filter))
      .flatMap(drawToPng(output))
      // For some reason the app does not exit without shutting down the EC
      .map(_ => ec.shutdown())
      .as(ExitCode.Success)
  }

  // Make a series using the datapoints in "file" that match "filter"
  def makeXYSeries(filter: String)(file: String)(implicit ec: ExecutionContext): IO[DrawMetrics.XYSeries] = {
    logger.info(s"Reading metric file at $file. Filtering for $filter")
    val path = Paths.get(file)

    io.file.readAll[IO](path, ec, 4096)
      .through(text.utf8Decode)
      .through(text.lines)
      .filter(_.contains(filter))
      .map(extractValue)
      .map(scaleY(filter))
      .zipWithIndex
      .map(_.swap)
      .compile
      .fold(new XYSeries(path.getFileName))({
        case (series, (x, y)) =>
          series.add(x.toDouble, y)
          series
      })
  }
  
  // Simply re-scale the values to MB if the filter is a memory filter
  def scaleY(filter: String)(value: Double) = filter match {
    case memory if memory.contains("memory") => value / MB
    case _ => value
  }

  // extract the value from a statsd packet to a Double
  def extractValue(line: String): Double = {
    val start = line.indexOf(':')
    val end = line.indexOf('|')
    line.substring(start + 1, end).toDouble
  }

  // Hide the X axis ticks which are meaningless here and set the title to the filter
  def configureChart(filter: String)(chart: XYChart) = {
    chart.plot.getDomainAxis.setTickLabelsVisible(false)
    chart.title_=(filter)
    chart
  }

  // Draw the chart to a file
  def drawToPng(output: String)(chart: XYChart) = IO {
    logger.info(s"Exporting graph to $output")
    chart.saveAsPNG(output)
    logger.info("Done")
  }
}
