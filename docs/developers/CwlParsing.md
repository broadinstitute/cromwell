# How CWL is Parsed

The ultimate entry point to parsing CWL is in `CwlV1_0LanguageFactory.validateNamespace`. 

Within that method, the process of parsing CWL can be thought of as performing these 4 steps:

1. Produce a canonical one-file representation of the CWL as a JSON object in memory
    - For the entrypoint into this process, see `CwlPreprocessor.preProcessCwl`
2. Use Circe to turn that representation into scala case classes
3. Convert the case classes into an appropriate WOM representations
4. Mix in inputs to produce a WOM Executable.

## Getting a canonical flat file

### From a source file and dependencies zip:

Note that the recursive implementation of `SALAD and flatten` in code can be found in: `CwlPreprocessor.preProcessCwl`.

The process looks like:

1. Write the workflow source file to disk
2. Make sure the dependencies are unzipped next to it
3. `SALAD and flatten` (see `CwlCanonicalizer.getCanonicalCwl`):
    - Run `cwltool --print-pre` against the source file to get a canonicalize JSON representation (aka `SALAD`) it.
        - NB: `cwltool --print-pre` must be able to resolve dependencies but it does not flatten them into the file.
    - Recursively `SALAD and flatten` any references into their own JSON representations of the contents.
    - Replace every reference in-place with the expanded JSON content.
4. At the end of this process we will have:
    - A single flat JSON 
    - Every imported step expanded to contain an in-line description
    - No more external references - everything required is provided in-line

### From a workflow URL:

The process looks similar to "source file and dependencies" except we never write anything locally:

1. `SALAD and flatten`:
    - Run `cwltool --print-pre` against the URL to get a canonicalize JSON representation (aka `SALAD`) it.
        - NB: `cwltool --print-pre` must be able to resolve dependencies but it does not flatten them into the file.
    - Recursively `SALAD and flatten` any references into their own JSON representations of the contents.
    - Replace every reference in-place with the expanded JSON content.
2. At the end of this process we will have:
    - A single flat JSON 
    - Every imported step expanded to contain an in-line description
    - No more external references - everything required is provided in-line


## Using Circe to get case classes

This process is more or less automatic - the case classes define the fields which they anticipate existing in the JSON.
and Circe does the rest.
    - Optional fields in case classes are allowed to not exist in the JSON
    - All fields in the JSON must be represented by fields in the appropriate case classes
    - In cases where the JSON may be structured in one of several ways, we use Shapeless coproducts to specify "one of" many options.

If a file is failing to parse the first stop should always be to double-check that the case classes accurately represent
the range of possible CWL JSON.

## Convert the case classes into WOM representations

- This is a recursive process of examining the CWL case classes and creating WOM values instead.
- Note that the CWL Language Factory (at least as of November 2018) does not implement the provision of WOM Bundles for 
imports. Instead, CWL Language Factory implements only the all-in-one method of converting CWL source and inputs files 
together into a WOM Executable.
    - Pragmatically, this just means that using the CWL language factory to process an import statement in a WDL file 
    is currently impossible (ie we cannot import CWL from WDL... yet).
