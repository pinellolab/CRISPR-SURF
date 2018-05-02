import os
import subprocess as sb
import sys

if len(sys.argv)==1:

	print	'\n- CRISPR-SURF Docker Container -\n\t- use `SURF` to use the command line\n\t- use `SURF_webapp` to start the web application\n'
	sys.exit(1)

if sys.argv[1]=='SURF':
	sb.call(["/opt/conda/bin/python", "/SURF/SURF.py"]+ sys.argv[2:])
elif sys.argv[1]=='SURF_webapp':
	sb.call(["/bin/bash", "-c", "/SURF/start_server_docker.sh"])
else:
	print	'\n- CRISPR-SURF Docker Container -\n\t- use `SURF` to use the command line\n\t- use `SURF_webapp` to start the web application\n'
