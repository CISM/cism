#!/usr/bin/awk -f
{if ($1 =="gl") print $11 " " ($12)}
