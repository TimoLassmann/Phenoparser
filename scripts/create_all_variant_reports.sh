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
usage: $0  -g <gemini_database> -d <phenoparser database> -t <report template>
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


          while getopts g:d:t:  opt
          do
              case ${opt} in
                  g) GEMINI_DATABASEPATH=${OPTARG};;
                  d) DATABASEPATH=${OPTARG};;
                  t) TEMPLATE=${OPTARG};;
                  *) usage;;
              esac
          done


          if [ "${GEMINI_DATABASEPATH}" = "" ]; then usage; fi
          if [ "${DATABASEPATH}" = "" ]; then usage; fi
          if [ "${TEMPLATE}" = "" ]; then usage; fi

          samples_array=$(gemini query -q "select name from samples" $GEMINI_DATABASEPATH)


          for i in ${samples_array[@]}; do
              step "Creating report for $i"
              try create_variant_report.sh  -i $i    -g $GEMINI_DATABASEPATH -d $DATABASEPATH  -t $TEMPLATE 
              next 
          done

          echo "DONE!!! Hurrah ";
      }

      main "$@"
