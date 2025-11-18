Per the spec,

"Provides input-specific hints in the form of an object. Each key within this
hint should refer to an actual input defined for the current task. A key may
also refer to a specific member of a struct/object input."

Currently, `localizationOptional` is the only hint Cromwell supports.

https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md#reserved-runtime-hints
https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md#localizationoptional
https://github.com/openwdl/wdl/blob/wdl-1.1/SPEC.md#inputs
