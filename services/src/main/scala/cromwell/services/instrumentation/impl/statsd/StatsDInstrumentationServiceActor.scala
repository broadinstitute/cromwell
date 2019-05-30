package cromwell.services.instrumentation.impl.statsd

import java.util.concurrent.{ConcurrentHashMap, TimeUnit}

import akka.actor.{Actor, ActorRef, Props}
import com.readytalk.metrics.{CromwellStatsD, StatsDReporter}
import com.typesafe.config.Config
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.instrumentation._
import cromwell.services.instrumentation.impl.statsd.StatsDInstrumentationServiceActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import nl.grons.metrics.scala.{DefaultInstrumented, Meter, MetricName}
import net.ceedubs.ficus.Ficus._

import scala.collection.JavaConverters._
import scala.concurrent.duration._

object StatsDInstrumentationServiceActor {
  val CromwellMetricPrefix: String = "cromwell"

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new StatsDInstrumentationServiceActor(serviceConfig, globalConfig, serviceRegistryActor))

  implicit class CromwellBucketEnhanced(val cromwellBucket: CromwellBucket) extends AnyVal {
    /**
      * Transforms a CromwellBucket to a StatsD path, optionally inserting a value between prefix and path 
      */
    def toStatsDString(insert: Option[String] = None) = (cromwellBucket.prefix ++ insert ++ cromwellBucket.path.toList).mkString(".")
  }
}

/**
  * Implementation of the Instrumentation service that sends metrics to a StatsD server.
  * Adapted from Workbench stack for StatsD instrumentation
  * Uses metrics-scala (https://github.com/erikvanoosten/metrics-scala) in combination with metrics-statsd (https://github.com/ReadyTalk/metrics-statsd)
  * metrics-scala is just a scala wrapper around dropwizard metrics which is a Java metric aggregation library (https://github.com/dropwizard/metrics)
  *
  * Metrics are stored and aggregated at the application level and periodically sent to a StatsD server via UDP.
  * Although StatsD supports a few different metric types, this infrastructure only send gauges that are pre-computed in a metrics registry.
  * This is not ideal because StatsD does not compute any statistics for gauges and simply forwards them to a backend.
  * It also adds some overhead on the application side that needs to compute its own statistics.
  * However it's simple working solution that is good enough for the current use case.
  *
  * If performance or statistics accuracy becomes a problem one might implement a more efficient solution
  * by making use of downsampling and / or multi metrics packets: https://github.com/etsy/statsd/blob/master/docs/metric_types.md
  */
class StatsDInstrumentationServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with DefaultInstrumented {
  val statsDConfig = StatsDConfig(serviceConfig)
  val cromwellInstanceIdOption: Option[String] = globalConfig.getAs[String]("system.cromwell_id")

  override lazy val metricBaseName = MetricName(CromwellMetricPrefix + cromwellInstanceIdOption.fold("")("." + _))

  val gaugeFunctions = new ConcurrentHashMap[CromwellBucket, Long]()

  StatsDReporter
    .forRegistry(metricRegistry)
    .prefixedWith(statsDConfig.prefix.orNull)
    .convertRatesTo(TimeUnit.SECONDS)
    .convertDurationsTo(TimeUnit.MILLISECONDS)
    .build(CromwellStatsD(statsDConfig.hostname, statsDConfig.port))
    .start(statsDConfig.flushRate.toMillis, TimeUnit.MILLISECONDS)

  override def receive = {
    case InstrumentationServiceMessage(cromwellMetric) => cromwellMetric match {
      case CromwellIncrement(bucket) => increment(bucket)
      case CromwellCount(bucket, value, _) => updateCounter(bucket, value)
      case CromwellGauge(bucket, value) => updateGauge(bucket, value)
      case CromwellTiming(bucket, value, _) => updateTiming(bucket, value)
    }
    case ShutdownCommand => context stop self
  }

  /**
    * We use a meter instead of a counter so we can at least get a rate (events / second)
    * The count is never reset though so the number keeps growing indefinitely until Cromwell is stopped / restarted
    */
  private def meterFor(bucket: CromwellBucket): Meter = {
    val name = bucket.toStatsDString()
    val counterName = metricBaseName.append(name).name
    metricRegistry.getMeters.asScala.get(counterName) match {
      // Make a new one if none is found
      case None => metrics.meter(name)
      // Otherwise use the existing one (after wrapping it in a metrics-scala object)
      case Some(meter) => new Meter(meter)
    }
  }

  /**
    * Increment the counter value for this bucket
    */
  private def increment(bucket: CromwellBucket) = meterFor(bucket).mark(1L)

  /**
    * Update the counter value for this bucket by adding (or subtracting) value
    */
  private def updateCounter(bucket: CromwellBucket, value: Long) = meterFor(bucket).mark(value)

  /**
    * Update the gauge value for this bucket
    */
  private def updateGauge(bucket: CromwellBucket, value: Long): Unit = {
    val newGauge = !gaugeFunctions.containsKey(bucket)
    gaugeFunctions.put(bucket, value)
    if (newGauge) {
      metrics.gauge(bucket.toStatsDString()){ gaugeFunctions.get(bucket) }
    }
    ()
  }

  /**
    * Adds a new timing value for this bucket
    */
  private def updateTiming(bucket: CromwellBucket, value: FiniteDuration) = {
    metrics.timer(bucket.toStatsDString()).update(value)
  }
}
