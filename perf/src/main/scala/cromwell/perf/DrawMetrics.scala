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
  val MB = 1024 * 1024
  
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
  
  def scaleY(filter: String)(value: Double) = filter match {
    case memory if memory.contains("memory") => value / MB
    case _ => value
  }

  def extractValue(line: String): Double = {
    val start = line.indexOf(':')
    val end = line.indexOf('|')
    line.substring(start + 1, end).toDouble
  }

  def configureChart(filter: String)(chart: XYChart) = {
    chart.plot.getDomainAxis.setTickLabelsVisible(false)
    chart.title_=(filter)
    chart
  }

  def drawToPng(output: String)(chart: XYChart) = IO {
    logger.info(s"Exporting graph to $output")
    chart.saveAsPNG(output)
    logger.info("Done")
  }
}
