# RZ Lab vs. SNARL Bd qPCR comparison README.md

## General instructions
Read README.md and and associated files before writing any code.
Never commit credentials to git, and warn user loudly if a given request appears to ask for this or may result in this.
Always check explicitly with user before committing large files (>5Mb) to git.
Always summarize proposed changes for user and seek confirmation before implementing
After implementing a cohesive set of changes, prompt user and assist in committing these changes to git/github and update README.md accordingly.

## Project-specific instructions
Default behavior is that Claude agents help write the R code, and user runs models/code. Agents can run small, non-intensive snippets of code if warented to debug, but user can generally do this and report back.