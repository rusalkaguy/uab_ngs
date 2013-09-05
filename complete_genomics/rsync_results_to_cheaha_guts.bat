@ECHO ""
@ECHO ##==== DS magic to get volume name - only from prompt! ====###
@ECHO vol ^| findstr "Volume.in.drive" ^> VOL.BAT
@ECHO for /f "usebackq delims=" %%i in (`type vol.BAT`) do set vol_fullstring=%%%%i
@ECHO set VOL_NAME=%%vol_fullstring:~22,40%%
@ECHO ECHO SET VOL_NAME=%%VOL_NAME%% ^> VOL.BAT
@ECHO ""
@ECHO ""
@ECHO push the data with cygwin magic
call vol.bat
SET TARGET=cheaha:/scratch/user/curtish/kimberly/cg/%VOL_NAME%
rsync -hmav "--include=**/" "--include=GS*-ASM.status.txt" "--include=*/*/*/ASM/*tsv*" "--include=*/*/*/ASM/REPORTS/*" "--exclude=*" . %TARGET%




