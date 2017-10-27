# Language Support

Find information on the workflow description languages that Cromwell supports here, as well as exciting hints about the future!


## Current Language Support

### WDL Draft 2
Cromwell started life as a WDL engine and that is our prefered execution language. For many examples on how to use WDL and some great getting-started resources you can see [the OpenWDL site](https://github.com/openwdl/wdl#getting-started-with-wdl).

Cromwell supports the majority of [Draft-2 of the WDL Spec](https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md).

#### Draft 2 Caveats

- Be careful when using `Object`. Cromwell has only half-hearted support and they are likely to be removed in future versions.
- Cromwell does not support nested `scatter`s.

## Future Language Support

### WDL Draft 3

We will be continuously building in support for [Draft-3](https://github.com/openwdl/wdl/tree/master/versions/draft-3) as the spec develops. 
Draft 2 support will be maintained during this process.

### CWL 1.0

From Cromwell v30 onwards, Cromwell will provide support for the CWL language, beginning with the core spec, and providing support for the following requirements:
- ShellCommandRequirement
- InlineJavascriptRequirement