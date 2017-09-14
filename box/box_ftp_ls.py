#!/usr/bin/env python

# see https://docs.python.org/2/library/ftplib.html#ftplib.FTP_TLS

import sys

from netrc import netrc, NetrcParseError
from getpass import getpass, getuser
from ftplib import FTP_TLS

def get_credentials(host):
    """
    Returns the ftp credentials for a host as a tuple containing 3 items in the order
        ('host', 'user', 'pass')
    This assumes a .netrc file exists in the users home directory and a valid host exists
    else you'll get prompted
    """

    try:
        infos = netrc().authenticators(host)
        # netrc returns (user,None,pass) - fix to match input for FTP
        credentials = (host,infos[0],infos[2])
    except NetrcParseError:
        print('Sorry no netrc for this host, please enter in your credentials')
        credentials = (host, getuser(), getpass())

    return credentials

def parse_line(line):
    """
    Takes a line of MLSD output, parses, and prints
    """
    d = dict((i.split("=")[0],i.split("=")[1]) for i in line.split(";")[0:4])
    d['Name']=line.split(";")[-1].lstrip()
    print(d)

def ftps_connect(credentials):
    """
    Takes in a tuple containing host, username and passwd and connects to ftp server and does stuff
    """
    cmd="LIST" # retrieves a list of files and information about those files.
    cmd="MLSD" # retrieves a machine readable list of files and information about those files

    # Using  *credentials unpacks your tuple for you WOWWW
    print('Connect to', credentials[0:2])
    try:
        ftps=FTP_TLS(*credentials)
        ftps.prot_p() # switch toi secure data connection
        # Do stuff here
        print(ftps.retrlines("MLSD"),parse_line)
        ftps.retrlines("MLSD",parse_line)
        #.....
    except:
        print('Sorry - FTPS connect failed:', credentials[0:1], "error:", sys.exc_info()[0])
        
    return 0



if __name__ == '__main__':
    # Always do input validation okay, please?

    # Size=497545189951;Modify=20170727181808.000;Create=20170621145109.000;Type=dir; Travis_Ptacek_Documentation
    hostname="ftp.box.com"
    credentials = get_credentials(hostname)
    rc=ftps_connect(credentials)
    sys.exit(rc)
