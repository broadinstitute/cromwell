# Language Support

Below are the Domain Specific Languages (DSL) that Cromwell currently supports and will soon support for describing your workflow.

## Current Language Support

### WDL Draft 2
Cromwell started life as a WDL engine and WDL draft-2 was our first language.
For many examples on how to use WDL and some great getting-started resources you can view [the OpenWDL site](https://github.com/openwdl/wdl#getting-started-with-wdl).

Cromwell supports the majority of [Draft-2 of the WDL Spec](https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md).

> *Known Issues:*
>
> - Be careful when using `Object`. They are superceded by 'struct' in WDL 1.0 and are being removed outright in WDL 2.0.
> - Cromwell does not support nested `scatter`s in draft-2.


### WDL 1.0

Cromwell supports WDL 1.0 and this is our most widely used version today.

As well as the changes to the WDL spec between draft-2 and 1.0, Cromwell also supports nested scatters and the [localization_optional](optimizations/FileLocalization.md) optimization in WDL 1.0.  

### WDL 1.1

Supported since version 92 in January 2026. Check out [the spec](https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md) to learn more.

### WDL 1.2

Supported since version 93. Cromwell implements [WDL 1.2.1](https://github.com/openwdl/wdl/releases/tag/v1.2.1). Check out [the spec](https://github.com/openwdl/wdl/blob/wdl-1.2/SPEC.md) to learn more.

The deprecated `runtime` section is still accepted, but it cannot be combined with the new `requirements` and `hints` sections in the same task.

> *Known Issues:*
>
> - `requirements`: `fpga` is not supported by any backend, and the WDL 1.2 `disks` formats (e.g. `"10 GiB"`, mount points) are passed to the backends unchanged. `container: "*"` is only recognized as a literal.
> - `hints`: `max_cpu`, `max_memory` and `short_task` are accepted but not acted on. The `disks`, `gpu` and `fpga` hints are ignored. Other hints are passed to the backend as runtime attributes, so Cromwell-specific attributes such as `preemptible` can be set as hints.
> - The `task` variable provides `name`, `id`, `container`, `cpu`, `end_time`, `meta`, `parameter_meta` and `ext`. Using `attempt`, `return_code`, `memory`, `gpu`, `fpga` or `disks` is a validation error.
> - `Int` values are 32-bit in Cromwell, not 64-bit as WDL 1.2 requires.
> - `size()` of a `Directory` is not supported yet.
> - Requirements and hints are overridden in the inputs JSON with `<task>.runtime` rather than `<task>.requirements` and `<task>.hints`.
> - WDL 1.2 is not supported on the PAPI v2 (Google Pipelines API) backend. Use GCP Batch instead.

## Future Language Support

### WDL 'development'

As the SPEC is being improved and honed, Cromwell continues to support the current `development` version of WDL. That
means that when (or shortly after) new versions are published, Cromwell will be ready to support them.
