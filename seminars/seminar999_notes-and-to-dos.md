Notes from R Markdown and Git(Hub) seminar
========================================================

notes to us for stuff to package and post ...

Prerequisite: install `knitr` (and its dependencies)

Flow problem: logically we should talk about Markdown before R Markdown but RStudio makes it easier to skip straight to R Markdown. Finessing this well is critical for making the steps clear, since various bits of the workflow differ for the two scenarios. For example, "Preview" vs. "Knit HTML", needing to load `markdown` vs `knitr` package, use of `markdownToHTML()` vs. `knit2HTML()`, etc.

Get a Markdown file, one (or more) of these ways:

  * Use the automagic skeleton RStudio provides. Do this: RStudio File > New File > R Markdown. Kind of jumping ahead because this is actual __R__ markdown. You can remove the "R chunks" to strip this down to plain Markdown.
  * Download a local copy of this [toy document](https://raw.github.com/jennybc/2013-11_sfu/master/simple-markdown.md) Jenny has posted on github.
  * Find one on the internet.
  * Write one for yourself!
  
Make sure the Markdown document's file name ends in `.md`. This is the canonical file extension. __Don't get creative here.__ If it's R Mardown, use `.rmd`.
  
Open this Markdown file in RStudio. It will appear in your editor pane.

Whenever you have a Markdown or R Markdown file open in the editor, the upper "bar" of the editor will have special powers.

Click on the question mark `?` and a Markdown Quick Reference will appear in another pane, by default in the same pane as help, plots, etc.

When you're ready, compile this Markdown file to HTML.

  * Via RStudio: Click on "Preview" (if editing Markdown) or "Knit HTML" (if editing R Markdown) to see your file after compilation to HTML.
  * Via R (in the Console or inside a script)
    - *now I wish I was writing in R markdown*
    - *fill this bit in*
  * Via the shell or in a Makefil
    - *fill this bit in*
    


Good place to learn more about git and to practice:
<try.github.io>

