#!@PYTHON@
#
#  chmconvert - Convert chemical networks into .chm format
#
#  Copyright (c) 2006-2011 Sebastien Maret
# 
#  This file is part of Astrochem.
#
#  Astrochem is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  Astrochem is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Astrochem.  If not, see <http://www.gnu.org/licenses/>.

import sys
import getopt
from libastrochem import network

def usage():
    """
    Display usage.

    """

    print """Usage: chmconvert [option] [file]

Common options:
   -h, --help         Display this help
   -V, --version      Display chmconvert version information
   -o, --output       Write the edited network in a file

See chmconvert(1) man page for a complete list of commands and options.
Report bugs to <sebastien.maret@obs.ujf-grenoble.fr>."""

def version():
    """
    Display version number.

    """

    print "This is chmconvert, version %s" % VERSION
    print """Copyright (c) 2006-2011 Sebastien Maret

This is free software. You may redistribute copies of it under the terms
of the GNU General Public License. There is NO WARRANTY, to the extent
permitted by law."""

def main():
    """
    Main function for chmconvert

    """
    
    # Parse options and check arguments. Display help message if
    # an unknown option is given.
    try:
	opts, args = getopt.getopt(sys.argv[1:], "hVo:",
				   ["help", "version", "output"])
    except getopt.GetoptError:
	usage()
	sys.exit(1)

    output = None
    for opt, arg in opts:
	if opt in ("-h", "--help"):
	    usage()
	    sys.exit()
	if opt in ("-V", "--version"):
	    version()
	    sys.exit()
	if opt in ("-o", "--output"):
	    output = arg

    # Check that we have a least a filename
    if len(args) == 1:
        network_file = args[0]
    else:
        sys.stderr.write("chmconvert: no file to convert.\n")
        sys.stderr.write("Type \'chmconvert --help\' for more information.\n")
        sys.exit(1)
        
    # Guess the format from the file extension
    if len(network_file.rsplit('.', 1)) == 2:
        network_file_base, network_file_ext = network_file.rsplit('.', 1)
    else:
        sys.stderr.write("chmconvert: file has no extension.\n")
        sys.exit(1)  
    if network_file_ext in ["osu", "kida"]:
        format = network_file_ext
    else:
        sys.stderr.write("chmconvert: unknown network format \"%s\".\n" 
                         % network_file_ext)
        sys.exit(1)
		
    # Open input and output files. We set fileout to stdout if
    # the --output option was not set.
    try:
        filein = open(network_file)
    except:
        sys.stderr.write("chmconvert: can't open %s.\n" % network_file)
        sys.exit(1)
    if not output:
        fileout = sys.stdout
    else:		
        try:
            fileout = open(output, 'w')
        except:
            sys.stderr.write("chmconvert: can't open %s.\n" % output)
            sys.exit(1)
	    
    # Read the network file and convert it
    try:
        net = network(filein, format = format)
    except Exception, err:
        sys.stderr.write("chmconvert: %s.\n" % err)
        sys.exit(1)        
    filein.close()

    # Write the converted file:
    net.write(fileout)
    fileout.close()
    
if __name__ == "__main__":
    main()
