# Language Support

Below are the Domain Specific Languages (DSL) that Cromwell currently supports and will soon support for describing your workflow.

## Current Language Support

### WDL Draft 2
Cromwell started life as a WDL engine and that is our recommended execution language. For many examples on how to use WDL and some great getting-started resources you can view [the OpenWDL site](https://github.com/openwdl/wdl#getting-started-with-wdl).

Cromwell supports the majority of [Draft-2 of the WDL Spec](https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md).

> *Known Issues:*
>
> - Be careful when using `Object`. Cromwell has only half-hearted support and they are likely to be removed in future versions.
> - Cromwell does not support nested `scatter`s.

## Future Language Support

### WDL Draft 3

We will be continuously building support for [Draft-3](https://github.com/openwdl/wdl/tree/master/versions/draft-3) as the spec develops.  
Draft 2 support will be maintained during this process.

### CWL 1.0

From Cromwell 30 onwards, Cromwell will provide support for Common Workflow Language (CWL), beginning with the core spec, and providing support for the following requirements:

* `ShellCommandRequirement`
* `InlineJavascriptRequirement`
