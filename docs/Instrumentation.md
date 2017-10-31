**Overview**

Cromwell's instrumentation support can be useful to collect utilization data in long-running, high-volume
production environments. While this instrumentation support can be used in smaller environments it will still require setting up a
[StatsD](https://github.com/etsy/statsd) server outside of Cromwell and it's possible not enough data would be produced to be useful.  

**StatsD**

Cromwell collects metrics while running and sends them to an internal service. The default implementation of this service
ignores these metrics, but Cromwell includes an alternate implementation that can forward metrics to a
[StatsD](https://github.com/etsy/statsd) server.
To specify this implementation in your configuration file:

```bash
services.Instrumentation.class = "cromwell.services.instrumentation.impl.statsd.StatsDInstrumentationServiceActor"
```
Make sure to configure your StatsD service:

```bash
services.Instrumentation.config.statsd {
    hostname = "localhost" # Replace with your host
    port = 8125 # Replace with your port
    # prefix = "my_prefix" # All metrics will be prefixed by this value if present.
    flush-rate = 1 second # Rate at which metrics are sent to the StatsD server
  }
```

There is also an additional configuration value that can be set: 

```bash
# Rate at which Cromwell updates its gauge values (number of workflows running, queued, etc...)
system.instrumentation-rate = 5 seconds
```

**Metrics**

The current StatsD implementation uses metrics-statsd to report instrumentation values.
metrics-statsd reports all metrics with a gauge type.
This means all metrics will be under the gauge section. We might add or remove metrics in the future depending on need and usage.
These are the current high level categories:

* `backend`
* `rest-api`
* `job`
* `workflow`
* `io`
