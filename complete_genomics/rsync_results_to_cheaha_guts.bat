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

echo starting... | email -V -r vera.dpo.uab.edu -f %USERNAME%@genome-curtis-h.ad.uab.edu -n "rsync_results_to_cheaha_guts.bat" -s "%PWD% starting rsync %VOL_NAME% to cheaha" %USERNAME%@uab.edu

rsync -hmav "--include=**/" "--include=manifest.all.out" "--include=GS*-ASM.status.txt" "--include=*/*/*/ASM/*tsv*" "--include=*/*/*/ASM/REPORTS/*" "--exclude=*" . %TARGET% > rsync.log.txt

type rsync.log.txt | email -V -r vera.dpo.uab.edu -f %USERNAME%@genome-curtis-h.ad.uab.edu -n "rsync_results_to_cheaha_guts.bat" -s "%PWD% completed rsync %VOL_NAME% to cheaha" %USERNAME%@uab.edu



