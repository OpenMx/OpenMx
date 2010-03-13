#!/usr/bin/python
import re
import sys
print "require(OpenMx)"
p = re.compile("\.\. code-block:: r\n\n(.+?)\n\n\w", re.DOTALL)
matches = p.finditer(sys.stdin.read())
for match in matches:
	print match.group(1)
