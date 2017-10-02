
# sync2box

Python program to backup a directory full of files to box.com via FTPS. 
 * compute MD5 sums for each file and store in .md5/<filename>
 * uses split to break files > 10G into 10G chunks before upload
 * uploads other files directly
 * scans up directory tree to find .ics.json to get ICS_ID (project Manager ID) and the correct relative path in box.com (see example below)
 * assumes your box external password/account are stored in ~/.netrc

## .ics.json 

```json
{
    "_comment": "Informatics Consulting Service (ICS) project metadata file"
    ,"ics_id":"xxx"
    ,"boxPath" : "consults/PI_NAME/PROOJECT_NAME"
}
```

## .netrc 

```
# for box.com
machine ftp.box.com
login BLAZERID@uab.edu
password BOX_EXTERNAL_PASSWORD
```

