# .bashrc

if [ -r /data ]; then
	# cheaha RC specific
	export PATH=~/bin:~/uab_ngs/slurm:~/uab_ngs/linux_plus:~/uab_ngs/lts:${PATH}
	#export MODULEPATH=/share/apps/ngs-ccts/modulefiles:${MODULEPATH}	
	export R_LIBS=~/R/R_LIB.rc
	# EasyBuild modules in my home dir
	export MYMODULEPATH=$(find ~/.local/easybuild/modules -maxdepth 1 -type d \! -name "modules" | perl -pe "s/\n/:/g")

	# DVC tools (git lfs)
	#alias git="module load Singularity/2.4.1-GCC-5.4.0-2.26 && singularity exec --bind /data /share/apps/ngs-ccts/simg/dvctools-0.2.simg git"
	#module load dvctools/latest	   # exciting version
	#module load dvctools/recommended   # stable version
	# make sure I don't forget
	alias git='echo module load dvctools "# because you forgot"; module load dvctools; alias lfsdiff="command singularity exec --bind /local --bind /data --bind /share --bind /scratch \$DVCTOOLS_SIMG git difftool -y --extcmd='\''/usr/bin/diff -c0 '\'' "; echo "now rerun: " git '
        # hack to fix dvctools git 
	alias lfsdiff='git'

	# same for snakemake
	alias snake='echo module purge \; module load shared rc-base snakemakeslurm; module purge; module load shared rc-base snakemakeslurm ; alias snakes="snakemakeslurm --use-conda"; alias snake="snakemake --use-conda "; alias snakes="snakemakeslurm --use-conda "; echo "now rerun: snake "'
	alias snakes='snake'
	
	alias getwd='cd $(pwd -LP)'

	# project things
	# kimberly
	alias tpt='cd /scratch/tptacek/pacbio_assembly_compiled/'
	alias hpgc='cd /scratch/curtish/kimberly/hpgc/genomes/svprobes/'
	alias rgc='cd /data/project/kimberly-lab/regeneron'
	alias cpc='cd /scratch/curtish/kimberly/cpc'
	alias guppy='cd /scratch/curtish/kimberly/ont_call_fast5'
	alias apr='cd /scratch/curtish/kimberly/apr'
	alias cser='cd /data/project/kimberly-cser'
	alias cpc='cd /scratch/curtish/kimberly/cpc'
	alias ontseq='cd /data/project/kimberly-lab/ont_seq'
	alias ontpipe='cd /data/project/kimberly-lab/ontpipe'
	alias ontpipet='tmux a -t ontpipe'
	alias pb='cd /data/project/kimberly-cser/pacbio_raw/'
	alias pbt='tmux a -t pb_raw'
	alias ichip='cd /data/project/ccts/user/curtish/kimberly/ics889_genotype_regn_dups/illumina_product_files'
	alias mk1c='ssh minit@138.26.142.125'

	# leal
	alias leal='cd /data/project/leallab/ics1515_genexus'
	alias leals='cd /data/project/leallab/ics1523_ont_spike/ics1523-Snakemake-ONT-Spike-Pipeline'
	alias lealz='cd /data/project/leallab/genexus_zips'
	alias lealt='tmux a -t leal'
	# markert/ HSV1
	alias hsv='cd /data/project/ccts/curtish/markert/ics1696_hsv1_nanopore'
	# ICTV
	alias vmr='cd /data/project/ccts/curtish/ictv/ICTVtaxablast'
	# Niederweis
        alias niederweis='cd /data/project/ccts/client/niederweis'

        # ICTV linux servers                                                                                                                                                                
        #                                                                                                                                                                                   
	alias ssh-ictv-app="ssh -i ~/.ssh/ICTV-curtish.pem ubuntu@app.ictv.global"
	alias ssh-ictv-test="ssh -i ~/.ssh/ICTV-curtish.pem ubuntu@test.ictv.global"
	alias ssh-ictv-prod="ssh -i ~/.ssh/ICTV-curtish.pem ubuntu@ictv.global"
	#                                                                                                                                                                                   
	# ONT Sequencer                                                                                                                                                                     
	#                                                                                                                                                                                   
	alias ssh-mk1c="ssh minit@138.26.142.125"
	#
	# linux short cuts
	#
	alias ffind='find . -type d -path "*/.*" -prune -o -not -name ".*" -type f'
	alias slurm_exports='(echo -e "\n\n\n\n\n# $host"; printenv | grep SLURM | awk "{print \"export \"\$0}"; egrep "SLURM(_NPROCS.*|_MEM_PER_CPU.*)*" )| tee ~/.$HOSTNAME; echo "# "source ~/.$HOSTNAME;'
else
	# cheaha UABGRID specific
	export PATH=~/uab_ngs/sge:~/uab_ngs/linux_plus:${PATH}
	export MYMODULEPATH=~/.modulefiles/curtish
	export R_LIBS=~/R/R_LIB.cheaha

fi
export MODULEPATH=${MYMODULEPATH}:$MODULEPATH

# env history
# curtish debugging 20150911
export ENV_HIST=${ENV_HIST}:/home/curtish/.bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Galaxy
#module load galaxy/galaxy-code-developer
export GALAXY_API_KEY=d6970f40cd5d725f92cf41d8af15b968
#----------------------------------------------------------------------
# /share/apps not mounted on several machines I accss, which causes hangs
#----------------------------------------------------------------------

# NGS tools
#if [ -f /share/apps/ngs-ccts/latest/bashrc ]; then
#        . /share/apps/ngs-ccts/latest/bashrc
#fi

# vcf-tools
#export PATH=/share/apps/ngs-ccts/vcftools_0.1.10/bin:$PATH
#export PERLLIB=/share/apps/ngs-ccts/vcftools_0.1.10/perl:$PERLLIB

#----------------------------------------------------------------------

# make sure backspace works in emacs
#stty erase 

# add parallel
export PATH=/home/$USER/bin/parallel:$PATH

# add ~/ics/qstat to path
export PATH=/home/$USER/uab_ngs/sge:$PATH
export PATH=/home/$USER/uab_ngs/utils:$PATH
export PATH=/home/$USER/uab_ngs/linux_plus:$PATH
#export PERLLIB=/home/$USER/lib/perl5/site_perl/5.8.8/

# add ~/lib to PYTHON path
export PYTHONPATH=/home/curtish/lib/python:$PYTHONPATH

# some aliases
alias fqp="ssh -X 172.20.100.45"
alias fqp-face="ssh -X face@172.20.100.45"
alias i2b2-ds="ssh -X 172.20.0.114"
alias facevm-01="ssh -X 172.20.0.115"
alias ls='\ls --color=tty'
#alias ls='lfs find . -maxdepth 1 | sort '
alias llst='\ls -lst --color=tty'
alias grep='\grep --color'
alias emacsn='emacs -nw '
alias gits='git status -s -uno'
alias gitp='git pull'
alias loadconda='PKGVER=$(grep "^[^#]*/data/user/curtish/Anaconda3/" ~/.condarc | cut -d / -f 5-); echo $PKGVER; module load $PKGVER; PROF=$EBROOTANACONDA3/etc/profile.d/conda.sh; if [ -e $PROF ]; then . $PROF; fi; PROF=$ANACONDA_ROOT/etc/profile.d/conda.sh; if [ -e $PROF ]; then . $PROF; fi'

# R
alias rnono="module load R/R-3.1.2; R --no-save --no-restore"

# TCOFFEE
#export DIR_4_TCOFFEE="/share/apps/ngs-ccts/t-coffee/Version_11.00.8cbe486"
#export MAFFT_BINARIES="$DIR_4_TCOFFEE/plugins/linux/"
#export CACHE_4_TCOFFEE="/home/curtish/.t_coffee/cache/"
#export TMP_4_TCOFFEE="$DIR_4_TCOFFEE/tmp/"
#export LOCKDIR_4_TCOFFEE="$DIR_4_TCOFFEE/lck/"
#export PERL5LIB="$PERL5LIB:$DIR_4_TCOFFEE/perl/lib/perl5"
#export EMAIL_4_TCOFFEE="chendrickson@uab.edu"
#export PATH="$DIR_4_TCOFFEE/bin:$PATH"

# personal PERL
export PERL_LOCAL_LIB_ROOT="$PERL_LOCAL_LIB_ROOT:/home/curtish/perl5";
export PERL_MB_OPT="--install_base /home/curtish/perl5";
export PERL_MM_OPT="INSTALL_BASE=/home/curtish/perl5";
export PERL5LIB="/home/curtish/perl5/lib/perl5:$PERL5LIB";
export PATH="/home/curtish/perl5/bin:$PATH";

# GIT LFS & Reproducible research
#export GIT_LFS_HOME=/share/apps/ngs-ccts/git-lfs/git-lfs-2.3.4
#export PATH=${GIT_LFS_HOME}:$PATH
#export QWRAP_SUPPORT_SCRIPTS_HOME=/share/apps/ngs-ccts/QWRAP-support-scripts/production
#export PATH=${QWRAP_SUPPORT_SCRIPTS_HOME}:$PATH

# The next line updates PATH for the Google Cloud SDK.
if [ -f '/home/curtish/google-cloud-sdk/path.bash.inc' ]; then source '/home/curtish/google-cloud-sdk/path.bash.inc'; fi

# The next line enables shell command completion for gcloud.
if [ -f '/home/curtish/google-cloud-sdk/completion.bash.inc' ]; then source '/home/curtish/google-cloud-sdk/completion.bash.inc'; fi

# time of last command in prompt
# see https://jakemccrary.com/blog/2015/05/03/put-the-last-commands-run-time-in-your-bash-prompt/

function timer_start {
  timer=${timer:-$SECONDS}
}

function timer_stop {
  timer_show=$(($SECONDS - $timer))
  unset timer
}

trap 'timer_start' DEBUG

if [ "$PROMPT_COMMAND" == "" ]; then
  PROMPT_COMMAND="timer_stop"
else
  PROMPT_COMMAND="$PROMPT_COMMAND; timer_stop"
fi

PS1='[last: ${timer_show}s][\w]$ '


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data/rc/apps/rc/software/Mamba/22.9.0-3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/rc/apps/rc/software/Mamba/22.9.0-3/etc/profile.d/conda.sh" ]; then
        . "/data/rc/apps/rc/software/Mamba/22.9.0-3/etc/profile.d/conda.sh"
    else
        export PATH="/data/rc/apps/rc/software/Mamba/22.9.0-3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

PS1='[last: ${timer_show}s][\w]$ '
