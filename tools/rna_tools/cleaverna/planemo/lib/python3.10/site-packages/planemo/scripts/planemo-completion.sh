#!/bin/sh

if [[ -n ${ZSH_VERSION-} ]]; then
    autoload -U +X bashcompinit && bashcompinit
fi

eval "$(_PLANEMO_COMPLETE=source planemo)"
