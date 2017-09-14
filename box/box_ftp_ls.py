#!/usr/bin/env python
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

def ftps_connect(credentials):
    """
    Takes in a tuple containing host, username and passwd and connects to ftp server and does stuff
    """

    # Using  *credentials unpacks your tuple for you WOWWW
    print('Connect to', credentials[0:2])
    try:
        ftps=FTP_TLS(*credentials)
        ftps.prot_p() # switch toi secure data connection
        # Do stuff here
        print(ftps.retrlines('LIST'))
        #.....
    except:
        print('Sorry - FTPS connect failed:', credentials[0:1], "error:", sys.exc_info()[0])
        
    return 0

if __name__ == '__main__':
    # Always do input validation okay, please?

    hostname="ftp.box.com"
    if len(sys.argv) > 1:
        hostname=sys.argv[1]

    credentials = get_credentials(hostname)

    rc= ftps_connect(credentials)
    sys.exit(rc)
