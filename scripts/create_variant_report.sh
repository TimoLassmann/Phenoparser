#!/usr/bin/env bash

export PATH=$PATH:~/gemini/anaconda/bin
export PATH=$PATH:~/bin
. /etc/init.d/functions

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


function file_exists ()
{
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



      pwd=$(pwd)

      function usage()
      {
cat <<EOF
usage: $0  -i <patient> -g <gemini_database> -d <phenoparser database> -t <template>
EOF
          exit 1;
      }

      main() {
          export PATH=$PATH:~/gemini/anaconda/bin
          export PATH=$PATH:~/bin

          PATIENT_ID=
          GEMINI_DATABASEPATH=
          DATABASEPATH=
          TEMPLATE=


          while getopts i:g:d:t:  opt
          do
          case ${opt} in
          i) PATIENT_ID=${OPTARG};;
      g) GEMINI_DATABASEPATH=${OPTARG};;
      d) DATABASEPATH=${OPTARG};;
      t) TEMPLATE=${OPTARG};;
      *) usage;;
      esac
      done

      if [ "${PATIENT_ID}" = "" ]; then usage; fi
      if [ "${GEMINI_DATABASEPATH}" = "" ]; then usage; fi
      if [ "${DATABASEPATH}" = "" ]; then usage; fi

      timestamp=$(date +"%m%d%y")

      reportname="Report_"$PATIENT_ID"_"$timestamp".Rmd";

      PATHNAME=$(dirname $TEMPLATE);
      echo "$PATHNAME";
      cat $TEMPLATE \
      | sed -e "s=VARPATIENT_ID=$PATIENT_ID=g" \
      | sed -e "s=VARGEMINI_DATABASEPATH=$GEMINI_DATABASEPATH=g" \
      | sed -e "s=VAROMIM_DATABASEPATH=$DATABASEPATH=g" \
      | sed -e "s=VARPHENOLYZER_DATABASEPATH=$PHENOLYZER_DATABASEPATH=g" \
      | sed -e "s=VARPATHTOHPOOBO=$PATHNAME=g" \
      > $reportname

      R -e "rmarkdown::render('$reportname')"


      echo "$timestamp $reportname";
   


      echo "DONE."
      }

      main "$@"
