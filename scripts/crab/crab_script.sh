#!/bin/bash

cmsRun -j FrameworkJobReport.xml -p PSet.py
echo '
============================================================
Analysing
============================================================
'
root -b -q crab_start.cxx crab_makeHistoss.cxx
