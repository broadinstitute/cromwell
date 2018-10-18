package cromwell.perf

import java.nio.file.Paths
import java.util.concurrent.Executors

import cats.effect.{ExitCode, IO, IOApp}
import cats.implicits._
import com.typesafe.scalalogging.StrictLogging
import fs2.{io, text}
import scalax.chart.XYChart

import scala.concurrent.ExecutionContext

object DrawMetrics extends IOApp with scalax.chart.module.Charting with StrictLogging {
  def run(args: List[String]): IO[ExitCode] = {
    implicit val ec = ExecutionContext.fromExecutorService(Executors.newFixedThreadPool(5))
    implicit val contextShift = IO.contextShift(ec)

    val filter = args.head
    val output = args.drop(1).head
    val sources = args.drop(2)

    sources
      .parTraverse[IO, IO.Par, DrawMetrics.XYSeries](makeXYSeries(filter))
      .map(XYLineChart(_))
      .map(configureChart)
      .flatMap(drawToPng(output))
      // For some reason the app does not exit without shutting down the EC
      .map(_ => ec.shutdown())
      .as(ExitCode.Success)
  }

  def makeXYSeries(filter: String)(file: String)(implicit ec: ExecutionContext): IO[DrawMetrics.XYSeries] = {
    logger.info(s"Reading metric file at $file. Filtering for $filter")

    io.file.readAll[IO](Paths.get(file), ec, 4096)
      .through(text.utf8Decode)
      .through(text.lines)
      .filter(_.contains(filter))
      .map(extractValue)
      .zipWithIndex
      .map(_.swap)
      .compile
      .fold(new XYSeries(filter))({
        case (series, (x, y)) =>
          series.add(x.toDouble, y)
          series
      })
  }

  def extractValue(line: String): Double = {
    val start = line.indexOf(':')
    val end = line.indexOf('|')
    line.substring(start + 1, end).toDouble
  }

  def configureChart(chart: XYChart) = {
    chart.plot.getDomainAxis.setTickLabelsVisible(false)
    chart
  }

  def drawToPng(output: String)(chart: XYChart) = IO {
    logger.info(s"Exporting graph to $output")
    chart.saveAsPNG(output)
    logger.info("Done")
  }
}
