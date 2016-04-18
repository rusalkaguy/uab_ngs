#!/bin/bash
tail -20000 /share/apps/galaxy/galaxy-latest/paster.log | grep -A 2 "^serving on http" | grep 2016
