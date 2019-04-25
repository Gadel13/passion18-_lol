#!/usr/bin/env sh
echo_green(){
    echo -e "\e[1;42m$1\e[0m"
}

set -e




cd Task4
echo "\n\n\n\n";
make cpplint
echo "\n"

 echo_green 'All tests cpplint done!'