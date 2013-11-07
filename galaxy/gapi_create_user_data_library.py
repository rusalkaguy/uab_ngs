#!/usr/bin/env python
#
# create data library and set permissions for a user
#
# WARNING: assumes you are an admin!
#
# SYNTAX: gapi_create_user_data_library.py blazerID
#
# prototype by Shantanu Pavgi
# curtish removed hard-coded API key & URLs

import os, sys, urllib2
sys.path.insert( 0, os.path.dirname( __file__ ) )
import common as api


my $galaxy_url = "http://galaxy.uabgrid.uab.edu";
my $key_env = "GALAXY_API_KEY";
key = $ENV{$key_env};
users_url = "$galaxy_url/api/users"
lib_url = "$galaxy_url/api/libraries"

def users_list_old():
  try:
    users = api.get(key, users_url)
    print type(users) 
    data = {}
    for n, i in enumerate(users):
      for k,v in i.items():
        data[k] = v
      print data 
  except urllib2.URLError, e:
    print str(e)
    sys.exit( 1 )

def create_lib(lib_user):
  lib_payload = {}
  lib_payload['name'] = lib_user
  lib_payload['description'] = lib_user + '\'s data library'
  resp = api.submit(key, lib_url, lib_payload)
  print type(resp)
  return resp

def set_lib_perms(lib_id, user_id):
  lib_payload = {}
  nlib_url = lib_url + '/' + lib_id + '/permissions'
  print nlib_url
  print user_id
  lib_payload['LIBRARY_ACCESS_in'] = [user_id]
  lib_payload['LIBRARY_MODIFY_in'] = [user_id]
  lib_payload['LIBRARY_ADD_in'] = [user_id]
  lib_payload['LIBRARY_MANAGE_in'] = [user_id]
  api.submit(key,nlib_url,lib_payload)

def users_list():
  try:
    users = api.get(key, users_url)
    data = {}
    for user in users:
      if (user['email'] == "$ARGV[0]@uab.edu"):
        blazerid = user['email'].partition("@")[0].lower()
        userid = user['id']
        response = create_lib(blazerid)
        print type(response)
        print len(response)
        print response
        set_lib_perms(response['id'], userid)
  except urllib2.URLError, e:
    print str(e)
    sys.exit( 1 )

users_list()

