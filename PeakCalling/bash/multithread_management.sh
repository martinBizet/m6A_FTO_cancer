############################################################
# Descrition: functions to management threads
# Date: 2016/01
# Author: Romy CHEN-MIN-TAO
############################################################

# MultiThread
# ===========

declare -a pids
function waitPids() {
    echo "Waiting for pids: ${pids[@]}"
    while [ ${#pids[@]} -ne 0 ]; do
        #echo "Waiting for pids: ${pids[@]}"
        local range=$(eval echo {0..$((${#pids[@]}-1))})
        local i
        for i in $range; do
            if ! kill -0 ${pids[$i]} 2> /dev/null; then
                echo "Done -- ${pids[$i]}"
                unset pids[$i]
            fi
        done
        pids=("${pids[@]}") # Expunge nulls created by unset.
        sleep 1
    done
    echo "Done!"
}       
 
function addPid() {
    desc=$1
    pid=$2
    echo "$desc -- $pid"
    pids=(${pids[@]} $pid)
}

function contener_command(){
    echo "[LOG] command: ${@}"
    $@
}

function free_processors(){
    nb_proc=$( cat /proc/cpuinfo | grep ^processor | wc -l )
    load_average=$( w | head -n 1 | awk -F " " '{gsub(",$","",$(NF-2)); gsub(",",".",$(NF-2));print $(NF-2)}' )
    free_proc=$( echo "${nb_proc} - ${load_average}" |bc -l ) 
    echo $free_proc
}

function check_load(){
    user_name=$1
    if [[ -z $user_name ]]
    then
        user_name=$USER
    fi
    # get number of processors
    nb_proc=$( cat /proc/cpuinfo | grep ^processor | wc -l )
    # get load average for 1 min
    load_average=$( w | head -n 1 | awk -F " " '{gsub(",$","",$(NF-2)); gsub(",",".",$(NF-2));print $(NF-2)}' )
    # get number free ram  
    ram_free=$( free -g | awk '{if ( $1 ~ /Mem/ ){print $4}}' )
    # get buffer 
    buffer_free=$( free -g | awk '{if ( $2 ~ /buffers/ ){print $4}}' )
    # get free swap
    swap_free=$( free -g | awk '{if ( $1 ~ /Swap/ ){print $4}}' )
    # job by user
    jobs_name=( $( ps -u $user_name | awk '{print $4}') )
    jobs_id=( $( ps -u $user_name | awk '{print $1}') )
    # report result
    echo "nb_proc:${nb_proc}\tload_average:${load_average}\tram_free:${ram_free}\tbuffer_free:${buffer_free}\tswap_free:${swap_free}\tjobs_name:${jobs_name[@]}"
}
