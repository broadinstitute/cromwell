# Language Support

Below are the Domain Specific Languages (DSL) that Cromwell currently supports and will soon support for describing your workflow.

## Current Language Support

### WDL Draft 2
Cromwell started life as a WDL engine and WDL draft2 was our first language!
For many examples on how to use WDL and some great getting-started resources you can view [the OpenWDL site](https://github.com/openwdl/wdl#getting-started-with-wdl).

Cromwell supports the majority of [Draft-2 of the WDL Spec](https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md).

> *Known Issues:*
>
> - Be careful when using `Object`. They are superceded by 'struct' in WDL 1.0 and are being removed outright in WDL 2.0.
> - Cromwell does not support nested `scatter`s in draft-2.


### WDL 1.0

Cromwell also supports WDL version 1.0.

As well as the changes to the WDL spec between draft-2 and 1.0, Cromwell also supports nested scatters and the [localization_optional](optimizations/FileLocalization.md) optimization in WDL 1.0.  


### CWL 1.0

Cromwell provides support for Common Workflow Language (CWL), beginning with the core spec, and most heavily used requirements.
If you spot a CWL feature that Cromwell doesn't support, please notify us using an issue on our github page!


## Future Language Support

### WDL 'development'

As the SPEC is being improved and honed, Cromwell continues to support the current `development` version of WDL. That
means that when (or shortly after) new versions are published, Cromwell will be ready to support them.
