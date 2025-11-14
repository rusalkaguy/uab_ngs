# ======================================================================
#
# [ICTV] host-colored terminal windows
#
# ======================================================================

# ICTV linux servers
#
source ~/ictv-zsh/ictv-zshrc
#
# Argonne/BV-BRC
#
alias ssh-argonne="$ICTV_ZSH_SET_GRASS; ssh -J ac.curtish@logins.cels.anl.gov ac.curtish@homes.cels.anl.gov; $ICTV_ZSH_SET_BASIC"
alias ssh-locust="$ICTV_ZSH_SET_GRASS; ssh -J ac.curtish@logins.cels.anl.gov ac.curtish@locust.mcs.anl.gov; $ICTV_ZSH_SET_BASIC"
alias ssh-peach="$ICTV_ZSH_SET_MANPAGE; ssh -J ac.curtish@logins.cels.anl.gov ac.curtish@peach.mcs.anl.gov; $ICTV_ZSH_SET_BASIC"
alias scp-peach="scp -J ac.curtish@logins.cels.anl.gov ac.curtish@peach.mcs.anl.gov "

#
# Cheaha cluster head nodew
#
alias ssh-rc="      $ICTV_ZSH_SET_NOVEL;           ssh ${USER}@cheaha.rc.uab.edu;      $ICTV_ZSH_SET_BASIC"
#
# ONT Sequencer
#
alias ssh-mk1c="    $ICTV_ZSH_SET_GRASS;           ssh minit@138.26.142.125;           $ICTV_ZSH_SET_BASIC"
alias ssh-mk1cSA="  $ICTV_ZSH_SET_SILVERAEROGEL;   ssh minit@138.26.142.125;           $ICTV_ZSH_SET_BASIC"
#
# Terminal
#
alias reset-zsh="$ICTV_ZSH_SET_BASIC"
alias ssh-reset-zsh="$ICTV_ZSH_SET_BASIC"

#
# pretty-print json files
#
alias jsonview="pygmentize -l json"

# ======================================================================
#
# Other settings/aliases
#
# ======================================================================

#
# short cuts
#
alias ffind='find . -type d -path "*/.*" -prune -o -not -name ".*" -type f'
alias vmr='cd ~/Documents/ICTV/VMR/VMR_to_BlastDB'
alias xls='open -a "Microsoft Excel" '
alias xlsx='open -a "Microsoft Excel" '

alias msl='cd ~/Documents/ICTV/load_msl'
alias msl40='cd ~/Documents/ICTV/load_msl40'
alias rgc='cd ~/Documents/kimberly/regeneron/fixit/analysis/pVCF/'

#
# other machines
#
alias ssh-nate='ssh NateAirBook' # 192.168.1.248 on ChezBlake

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/homebrew/Caskroom/miniforge/base/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh" ]; then
        . "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh"
    else
        export PATH="/opt/homebrew/Caskroom/miniforge/base/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

export PATH=/Users/curtish/edirect:${PATH}
