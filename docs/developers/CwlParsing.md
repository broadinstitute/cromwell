# How CWL is Parsed

The process of parsing CWL runs approximately like this:

1. Produce a canonical one-file representation of the CWL as a JSON object in memory 
1. Use Circe to turn that representation into scala case classes
1. Convert the case classes into an appropriate WOM representations

## Getting a canonical flat file

### From a source file and dependencies zip

The process looks like:

1. Write the workflow source file to disk
1. Make sure the dependencies are unzipped next to it
1. SALAD and flatten:
    - Run `cwltool --print-pre` against the source file to get a canonicalize JSON representation (aka `SALAD`) it.
        - NB: `cwltool --print-pre` must be able to resolve dependencies but it does not flatten them into the file.
    - Recursively SALAD and flatten and references into their own JSON representations of the contents.
    - Replace every reference in-place with the expanded JSON content.
1. At the end of this process we will have:
    - A single flat JSON 
    - Every imported step expanded to contain an in-line description
    - Only relative names (ie we remove all references to `file://tmp/...` since those are not deterministic during Cromwell restarts)

## Using Circe to get case classes

This process is more or less automatic - the case classes define the fields which they anticipate existing in the JSON.
and Circe does the rest.
    - Optional fields in case classes are allowed to not exist in the JSON
    - All fields in the JSON must be represented by fields in the appropriate case classes
    - In cases where the JSON may be structured in one of several ways, we use Shapeless coproducts to specify "one of" many options.

If a file is failing to parse the first stop should always be to double-check that the case classes accurately represent
the range of possible CWL JSON.

## Convert the case classes into WOM representations

As of today, CWL doesn't produce WOM bundles. Instead, it goes straight to executable.
