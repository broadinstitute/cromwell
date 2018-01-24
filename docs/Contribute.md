**Contribution Guidelines**

Cromwell is an open-source project and we heartily welcome community contributions to both our code and our documentation. Here are some guidelines for adding documentation and recommendations on where we could use help the most. 

First off, here are useful links:

* Cromwell documentation: [cromwell.readthedocs.io/en/develop](http://cromwell.readthedocs.io/en/develop)
* Source on Github: [github.com/broadinstitute/cromwell](https://github.com/broadinstitute/cromwell/tree/develop/docs)
* Builds on ReadTheDocs: [readthedocs.org/projects/cromwell](https://readthedocs.org/projects/cromwell/builds/)
* Information on how to build and view the documentation locally: [mkdocs.readthedocs.io](https://mkdocs.readthedocs.io/en/stable/#installation)

### Formatting

The documentation is written in Markdown. Click here for a [Github Guide on Markdown](https://guides.github.com/features/mastering-markdown/). Click here for more [tips from MkDocs](http://www.mkdocs.org/user-guide/writing-your-docs/).

### Styling

**Links:**

* Relative path: \[link text](Destination_Page) without the .md extension  
	* Ex: \[Getting Started](GettingHelp) produces [Getting Started](GettingHelp) 
* Absolute path: \[link text](www.broadinstitute.com)  
	* Ex: \[Broad Institute](www.broadinstitute.org) produces [Broad Institute](www.broadinstitute.org)
* Anchor links: \[anchor text](../Path/To/Page#Anchor)  
	* Ex: \[HPC filesystems](backends/HPC#filesystems) produces [HPC filesystems](backends/HPC#filesystems)

**Left-side menu:**  
To add or remove items from the menu, edit [mkdocs.yml](https://github.com/broadinstitute/cromwell/blob/develop/mkdocs.yml) in Cromwell.

**Code:**  
To style a word of code, use \` before and after the word.  
For example: \`file.json\` to produce `file.json`. 

To style a block of code, use \`\`\` before and after the block of code.  
Ex:



To use syntax highlighting, include the language after the first three \`\`\`  
For example:   

	```wdl  
		workflow myWorkflow {
    call myTask
	}   
	```	
to produce  
```bash  
	workflow myWorkflow {
    call myTask
}   
```	

Images
Relative path: 
Absolute path: 

REST API:
Edit the cromwell.yaml to make any changes to the REST API.
The REST API markdown may be regenerated from the cromwell.yaml by running
sbt generateRestApiDocs
Commit the changes to both the cromwell.yaml and the generated RESTAPI.md



### Writing Tips

1. Keep it clear, accurate, and concise, 
2. Put the most important information first,
3. Use the second person, use “you” instead of “the user”,
4. No passive verbs (everything is done by something),
5. Link to the original source, don't repeat documentation. 

### FAQs

**_Why isn't my documentation showing up?_**

**Have you confirmed that your PR is merged?**  
If not, [go merge it](https://github.com/broadinstitute/cromwell/pulls). This will trigger an automatic build.  

**Has the build finished?**   
[Check build status here](https://readthedocs.org/projects/cromwell/builds/).

**Did you add the file(s) to the YAML file?**  
If not, [add it here](https://github.com/broadinstitute/cromwell/blob/develop/mkdocs.yml).