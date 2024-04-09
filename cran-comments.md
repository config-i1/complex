---
title: "Cran Comments"
author: "Ivan Svetunkov"
date: "9 April 2024"
output: html_document
---

# Resubmission comment
Fixed issues outlined by Konstanze Lauseker, including "dontrun" for a function and a monograph details in the "DESCRIPTION" file.


# Version
This is the initial release of the package `complex`, v1.0.0


# Test environments
* Local Ubuntu 23.10, R 4.3.1
* Github Actions
* win-builder (devel and release)
* rhub with rhub::check_for_cran() command


## Local Ubuntu: R CMD check results
>R CMD check results
>0 errors | 0 warnings | 0 notes


## Github actions
Successful checks for:

- Windows latest release (2022 10.0.20348) with R 4.3.3
- MacOS latest release (macOS 12.7.4) with R 4.3.3
- Ubuntu 22.04.4 with R 4.3.3


## win-builder
>Possibly misspelled words in DESCRIPTION:
>  Sergey (13:27)
>  Springer (14:35)
>  Svetunkov (13:17, 13:38)

They are not misspelled.


## R-hub
### Windows Server 2022, R-devel, 64 bit
>* checking for non-standard things in the check directory ... NOTE
>Found the following files/directories:
>  ''NULL''

This is absurd! There is no file/directory "NULL" in the project. This seems like an error on the server side.

>* checking for detritus in the temp directory ... NOTE
>Found the following files/directories:
>  'lastMiKTeXException'

Once again, no such file or directory present in the package, probably an error on the server.

# Downstream dependencies
No reverse dependencies yet
