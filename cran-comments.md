---
title: "Cran Comments"
author: "Ivan Svetunkov"
date: "8 April 2024"
output: html_document
---

## Version
This is the initial release of the package `complex`, v1.0.0


## Test environments
* local Ubuntu 23.10, R 4.3.1
* win-builder (devel and release)
* rhub with rhub::check_for_cran() command


## R CMD check results
>R CMD check results
>0 errors | 0 warnings | 0 notes


## Github actions
Successful checks for:

- Windows latest release with R 4.2.3
- MacOS latest macOS Big Sur 10.16 with R 4.2.3
- Ubuntu 20.04.5 with R 4.2.3


## win-builder



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

## Downstream dependencies
No reverse dependencies yet
