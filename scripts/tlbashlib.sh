#!/bin/bash

# thanks to http://tech.franzone.blog/2008/08/25/bash-script-init-style-status-message/
# Column number to place the status message
RES_COL=60
# Command to move out to the configured column number
MOVE_TO_COL="echo -en \\033[${RES_COL}G"
# Command to set the color to SUCCESS (Green)
SETCOLOR_SUCCESS="echo -en \\033[1;32m"
# Command to set the color to FAILED (Red)
SETCOLOR_FAILURE="echo -en \\033[1;31m"
# Command to set the color back to normal
SETCOLOR_NORMAL="echo -en \\033[0;39m"

# Function to print the SUCCESS status
echo_success() {
  $MOVE_TO_COL
  echo -n "["
  $SETCOLOR_SUCCESS
  echo -n $"  OK  "
  $SETCOLOR_NORMAL
  echo -n "]"
  echo -ne "\r"
  return 0
}

# Function to print the FAILED status message
echo_failure() {
  $MOVE_TO_COL
  echo -n "["
  $SETCOLOR_FAILURE
  echo -n $"FAILED"
  $SETCOLOR_NORMAL
  echo -n "]"
  echo -ne "\r"
  return 1
}

step() {
    echo -n "STEP: $@"

    STEP_OK=0
    [[ -w /tmp ]] && echo $STEP_OK > /tmp/step.$$
}

try() {
    # Check for `-b' argument to run command in the background.
    local BG=

    [[ $1 == -b ]] && { BG=1; shift; }
    [[ $1 == -- ]] && {       shift; }

    # Run the command.
    if [[ -z $BG ]]; then
        "$@"
    else
        "$@" &
    fi

    # Check if command failed and update $STEP_OK if so.
    local EXIT_CODE=$?

    if [[ $EXIT_CODE -ne 0 ]]; then
        STEP_OK=$EXIT_CODE
        [[ -w /tmp ]] && echo $STEP_OK > /tmp/step.$$

        if [[ -n $LOG_STEPS ]]; then
            local FILE=$(readlink -m "${BASH_SOURCE[1]}")
            local LINE=${BASH_LINENO[0]}

            echo "$FILE: line $LINE: Command \`$*' failed with exit code $EXIT_CODE." >> "$LOG_STEPS"
        fi
    fi

    return $EXIT_CODE
}

next() {
    [[ -f /tmp/step.$$ ]] && { STEP_OK=$(< /tmp/step.$$); rm -f /tmp/step.$$; }
    [[ $STEP_OK -eq 0 ]]  && echo_success || echo_failure
    echo

    return $STEP_OK
}


file_exists () {
    if ! [ "$1" ]
    then
        echo "file exists function needs an input";
        return 1;
    fi
    if [[ -f "$1" ]]; then
        return 0;
    else
        return 1;
    fi;
}
