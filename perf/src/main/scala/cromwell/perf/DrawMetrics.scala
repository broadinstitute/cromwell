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
    val blockingExecutionContext = ExecutionContext.fromExecutorService(Executors.newFixedThreadPool(2))
    val filter = args(2)

    logger.info(s"Reading metric file at ${args.head}. Filtering for $filter")

    io.file.readAll[IO](Paths.get(args.head), blockingExecutionContext, 4096)
      .through(text.utf8Decode)
      .through(text.lines)
      .filter(_.contains(filter))
      .map(extractValue)
      .zipWithIndex
      .map(_.swap)
      .compile
      .fold(new XYSeries(args(2)))({
        case (series, (x, y)) =>
          series.add(x.toDouble, y)
          series
      })
      .map(XYLineChart(_))
      .map(configureChart)
      .flatMap(drawToPng(args(1)))
      // For some reason the app does not exit without shutting down the EC
      .map(_ => blockingExecutionContext.shutdown())
      .as(ExitCode.Success)
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
