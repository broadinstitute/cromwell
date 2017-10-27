_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the Instrumentation page?  
*why do they care about instrumentation and metrics? what do they use it for?*
2. What do they need to know first?  

3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

**Overview**

Cromwell's instrumentation support can be useful to collect utilization data in long-running, high-volume
production environments. While this instrumentation support can be used in smaller environments it will still require setting up a
[StatsD](https://github.com/etsy/statsd) server outside of Cromwell and it's possible not enough data would be produced to be useful.  


**StatsD**

Cromwell collects metrics while running and sends them to an internal service. The default implementation of this service
ignores these metrics, but Cromwell includes an alternate implementation that can forward metrics to a
[StatsD](https://github.com/etsy/statsd) server.
To specify this implementation in your configuration file:

```text
services.Instrumentation.class = "cromwell.services.instrumentation.impl.statsd.StatsDInstrumentationServiceActor"
```
Make sure to configure your StatsD service:

```text
services.Instrumentation.config.statsd {
    hostname = "localhost" # Replace with your host
    port = 8125 # Replace with your port
    # prefix = "my_prefix" # All metrics will be prefixed by this value if present.
    flush-rate = 1 second # Rate at which metrics are sent to the StatsD server
  }
```

There is also an additional configuration value that can be set: 

```text
# Rate at which Cromwell updates its gauge values (number of workflows running, queued, etc...)
system.instrumentation-rate = 5 seconds
```

**Metrics**

The current StatsD implementation uses metrics-statsd to report instrumentation values.
metrics-statsd reports all metrics with a gauge type.
This means all metrics will be under the gauge section. We might add /remove metrics in the future depending on need and usage.
These are the current high level categories:

* `backend`
* `rest-api`
* `job`
* `workflow`
* `io`
