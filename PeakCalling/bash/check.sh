# Check
# =====

# Test if variables exist, if not quit program
function is_empty(){
    name_value=$1
    value=$2
    if [[ -z $value ]]
    then
        echo "[ERROR] variable ${name_value} is empty"
        exit 1
    fi
}

# Test if directory exist, if yes move this directory on name_of_directory + "_OLDER"+ global variable DATE
# ... and create this directory empty 
function dir_exist(){
    local directory=$1
    dir_is_full=$( ls ${directory} )
    echo $dir_is_full
    if [[ ! -d $directory ]]
    then
        echo "[LOG] create directory: ${directory}" >&2
        mkdir -p ${directory}
    elif [[ ! -z $dir_is_full ]]
    then
        echo "[LOG] create new directory: ${directory}_OLDER${DATE}" >&2
        mkdir -p ${directory}_OLDER${DATE}
        echo "[LOG] move the content of ${directory}/* on new directory: ${directory}_OLDER${DATE}" >&2
        mv ${directory}/* ${directory}_OLDER${DATE}
    fi
}

# Describe path
# ... return :
# "f" if is a regular file 
# "s" if is not zero size
# "r" if file has read permission
# "w" if file has write permission 
# "x" if file has execute permission 
# or return combinaison of letters if file is in several cases
function describe_file(){ 
    # arguments
    local fpath=$1
    local test_file=$2
    # variables
    local describe_file=""
    if [[ -f $fpath ]]
    then
        describe_file="${describe_file}f"
    else
        echo "[ERROR] file doesn't exist: ${fpath}" >&2
    fi
    if [[ -s $fpath ]]
    then
        describe_file="${describe_file}s"
    else
        echo "[WARNING] file is empty : ${fpath}" >&2
    fi
    if [[ -r $fpath ]]
    then
        describe_file="${describe_file}r" 
    else
        echo "[WARNING] not permission to read : ${fpath}" >&2
    fi
    if [[ -w $fpath ]]
    then                            
        describe_file="${describe_file}w" 
    else                            
        echo "[WARNING] not permission to write in file : ${fpath}" >&2
    fi 
    if [[ -x $fpath ]]
    then                            
        describe_file="${describe_file}x" 
    else                            
        echo "[WARNING] not permission to execute : ${fpath}" >&2
    fi
    return $describe_file
}

# Test if you have permission to run one exe and if it exist 
# ... if not stop program
function exe_exist(){
    local fexe=$1
    if [[ -z $fexe ]]
    then
        echo "[ERROR] miss argument" >&2
        exit 1
    elif [[ ! -x $fexe ]]
    then                            
        echo "[ERROR] not permission to execute : ${fexe}" >&2
        exit 1
    fi 
}

# Test if file exist and if not empty
# ... if not stop program
function f_exist(){
    local fpath=$1
    if [[ -z $fpath ]]             
    then                          
        echo "[ERROR] miss argument" >&2                                                                                                      
        exit 1                    
    elif [[ ! -s $fpath ]]
    then
        echo "[ERROR] file is empty : '${fpath}'" >&2
        exit 1
    fi
}

# If file exist and not empty 
# ... move this file on ${fpath}_OLDER${DATE} with DATE is global variable
# ... and create new file empty
function create_f(){
    fpath=$1
    if [[ -s $fpath ]]   
    then
        mv $fpath ${fpath}_OLDER${DATE}
    fi
    > $fpath
}

# Test a list of file from pattern
function fpattern_exists(){
    local fpattern=$1
    list_f=$( ls $fpattern )
    is_empty fpattern $fpattern
    for f in ${list_f[@]}
    do
        echo "[LOG] test file: ${f}" >&2
        f_exist "${f}"
    done
}

# Test if function exist
function fn_exists() {
    local fn=$1
    test_is=$( type -t $fn )
    if [[ "${test_is}" != 'function' ]]
    then
        echo "[ERROR] This function doesn't exist : ${fn}" >&2
        exit 1
    fi
}

# Test if last action is OK
function action_ok(){
    # run command
    echo "[LOG] command to check: ${@}" >&2
    ${@}
    test_ok=$?
    # test is ok
    if [[ $test_ok -eq 1 ]]
    then
        echo "[ERROR] last action" >&2
        exit 1
    fi
}

# Test if global variables exist in bash environnement:
function check_env(){    
    echo "[LOG] check environnement" >&2
    # Arguments
    local list_variable=( $@ )
    # Scan arguments
    for d in ${list_variable[@]}
    do
        eval "[[ ! -z \${$d} ]]"     
        result=$?
        if [[ $result -eq 1 ]]
        then
            echo "[ERROR] variable ${d} don't exist" >&2
            echo 1
            return 1
        fi
    done
    echo 0
    return 0
}
